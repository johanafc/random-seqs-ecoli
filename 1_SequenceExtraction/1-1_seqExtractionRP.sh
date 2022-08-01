#!/bin/bash
#  Pipeline for analysis of random peptide effect of e-coli in experiments
# described in Neme, 2017.
## AUTHOR: JFajardo (04.2017)
## EDITED: 06.2020

## USAGE: bash ./1-1_seqExtractionRP.sh ARG
# (ARG= Experiment number)

## Sequencing fastq files have the format A_BC.fastq
## A. 1-1 = day(or time)-replicate
## B. _S2_L001_ = Randomly assigned Illumina sequencing file name
## C. R1 or R2 = Forward or Reverse reads, respectively

function usage {
	echo "Usage: $(basename $0) [-fdta]" 2>&1
	echo '   -f  Experiment folder name'
	echo '   -d  Path to working directory with experiment folders and scripts'
	echo '   -t  Path to trimmomatic executable'
	echo '   -a  Trimmomatic adapter file'
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

while getopts ':fdta' OPTION; do
	case $OPTION in
		f)
			exp="${OPTARG}/"
			echo "Experiment folder name: $exp"
			;;
		d)
			main_path="${OPTARG}/"
			echo "Working directory: $main_path"
			;;
		t)
			trimmomatic="${OPTARG}"
			echo "Path to trimmomatic: $trimmomatic"
			;;
		a)
			path_adapters="${OPTARG}"
			echo "Using trimming adapters: $path_adapters"
			;;

		?)
      		echo "Invalid option: -${OPTARG}."
      		echo
      		usage
      ;;
	esac
done

path_exp=${main_path}${exp}"/"

cd ${path_exp}
printf "%b\n" "Starting analysis "${exp}
mkdir ${path_exp}results
mkdir ${path_exp}TempReads

##First step: Generate a name table where $1= A, $2=Full_fastq_file_name(R1),
# $3=Full_fastq_file_name(R2)
# Output (1 file): names_exp-#.tsv
printf "%b\n" "-----------------------------------------------------------------------------"
printf "%b\n" "Step 1: Generating name table for ${exp}"
printf "%b\n" "-----------------------------------------------------------------------------"
printf '%(%a,%D;%T)T.\n' -1

ls -v *.fastq.gz |
	grep R2 |
	awk '{print $1" "$1" "$1}' |
	sed 's/_/ /' |
	awk '{print $1" "$3" "$4}' |
	sed 's/R2/R1/' |
	awk '{print $1"\t"$2"\t"$3"\t"$2"\t"$3}'|
	sed -e 's/.gz//' -e 's/.gz//'>names_${exp}.tsv
printf "%b\n" "Name table generated for ${exp}"


# Second step: Generate a non redundant database of peptides for each experiment
# Output (4 x #ofExperiments): $Exp.fa; $Exp-pre.org.gz; $Exp.orf; $Exp.pep
no=0
cat names_${exp}.tsv|
	while read -r rep r1 r2 r1gz r2gz;
	do
		no=$((${no}+1))
		printf "%b\n" "-----------------------------------------------------------------------------"
		printf "%b\n" "Step 2-${no}: Getting ${exp}, ${rep} single dereplicated files for peptide database"
		printf "%b\n" "-----------------------------------------------------------------------------"
		printf '%(%a,%D;%T)T.\n' -1

		java -jar ${trimmomatic} PE -threads 24 -phred33 ${r1gz} ${r2gz} "PE-"${r1} "SE-"${r1} "PE-"${r2} "SE-"${r2} ILLUMINACLIP:${path_adapters}TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 AVGQUAL:20

## Merge paired reads using usearch (V10.0) and modify sequence names to prevent
## them from being shortened by getorf
		printf "%b\n" "Merging paired-end reads"
		usearch10 -fastq_mergepairs "PE-"${r1} -reverse "PE-"${r2} -fastaout ${rep}.tmp -threads 24 -fastq_maxdiffs 30 -fastq_minmergelen 100
		seqtk seq ${rep}.tmp |
		sed -e 's/^>..*-\(..*\) ..*$/>\1/g' -e 's/:/_/g' >${rep}.merged
		printf "%b\n" "-Done- Merging paired-end reads"

## Keep only reads that have the M13 primer and the FLAG features from the plasmid
		perl ${main_path}1-1_rmSTART.pl ${rep}.merged ${main_path}1-1_Primers.txt

## Get nucleotide sequences of all ORFs between a start and a stop codon. The
## common transformation of seqtk seq removes extra line breaks in the fasta files
		printf "%b\n" "-----------------------------------------------------------------------"
		printf "%b\n" "Step 3-${no}: Getting random peptide ORFs from reads ${exp}, rep ${rep}"
		printf "%b\n" "-----------------------------------------------------------------------"
		printf '%(%a,%D;%T)T.\n' -1

		getorf -sequence ${rep}.merged.clean.fa -outseq ${rep}.pre.orf -minsize 12 -find 3
		sed 's/ .*//g' ${rep}.pre.orf>tmp

		seqtk seq tmp |
		grep -A1 "_1"$ |
		sed "/--/d" >${rep}.orf

		rm ${rep}.pre.orf tmp
		printf "%b\n" "
		ORFs saved to ${rep}.orf in fasta format"

		printf "%b\n" "-----------------------------------------------------------------------"
		printf "%b\n" "Step 4-${no}: Dereplicating full length clean reads and individual ORF sequences"
		printf "%b\n" "-----------------------------------------------------------------------"
		printf '%(%a,%D;%T)T.\n' -1

		##De-replicate full, clean reads
		usearch10 -fastx_uniques ${rep}.merged.clean.fa -sizeout -threads 24 -minuniquesize 2 -fastaout ${rep}".merged.clean.drp.fa"
		printf "%b\n" "Temporal file "${rep}".merged.clean.drp.fa"

		##De-replicate ORFS
		usearch10 -fastx_uniques ${rep}.orf -sizeout -minuniquesize 2 -threads 24 -fastaout ${rep}"_drp.orf"
		printf "%b\n" "Temporal file "${rep}"_drp.orf created"
	done

## Delete and compress intermediate files
mv SE* TempReads/.
mv PE* TempReads/.
rm *.tmp
gzip *.fastq

## Generate a dereplicated file of all cycles and replicates using USEARCH10:
## The 32bit version of USEARCH10 only allows to use up to 4Gb of memory
## at the time, if one tries to dereplicate a bigger file, it fails. To prevent
## this, we have to split the files in 3K lines, dereplicate each split file and
## merge.

#Full seqs
cat *.merged.clean.drp.fa > all.drp_clean.fa
faToTab -type=protein all.drp_clean.fa tmp
wc -l tmp
sort -k 2 -d tmp |
	awk '{print ">"$1"\n"$2}'>tmp.fa
split -d -l 3000000 tmp.fa split_
rm tmp*
for spl1 in split_*;
do
	usearch10 -fastx_uniques ${spl1} -sizein -sizeout -minuniquesize 4 -threads 24 -fastaout "drp_"${spl1}
done
for i in `seq 1 3`;
do
	cat drp_split_* >drp${i}_all.fa
	rm drp_split_*
	wc -l drp${i}_all.fa
	faToTab -type=protein drp${i}_all.fa tmp
	wc -l tmp
	sort -k 2 -d tmp |
		awk '{print ">"$1"\n"$2}'>drp${i}_all.fa
	split -d -l 3000000 drp${i}_all.fa split_
		for spl1 in split_*;
		do
			usearch10 -fastx_uniques ${spl1} -sizein -sizeout -minuniquesize 4 -threads 24 -fastaout "drp_"${spl1}
		done
	rm split_* tmp*
done

cat drp_split_* > ${exp}".clean_drp.fa"
rm all.drp_clean.fa

printf "%b\n" "-----------------------------------------------------------------------"
printf "%b\n" "Dereplicated full sequence files generated:"
wc -l ${exp}"_drp.clean.fa"
printf "%b\n" "-----------------------------------------------------------------------"
printf '%(%a,%D;%T)T.\n' -1

printf "%b\n" "-----------------------------------------------------------------------"
printf "%b\n" "Step 5: Generating dereplicated file of ORFs and full length reads in "${exp}
printf "%b\n" "-----------------------------------------------------------------------"
printf '%(%a,%D;%T)T.\n' -1

cat *_drp.orf > all_drp_orf.fa
faToTab -type=protein all_drp_orf.fa tmp
sort -k 2 -d tmp |
	awk '{print ">"$1"\n"$2}'>tmp.fa
split -d -l 3000000 tmp.fa split_orf_
rm tmp*
for spl1 in split_orf_*;
do
	usearch10 -fastx_uniques ${spl1} -sizein -sizeout -minuniquesize 4 -threads 24 -fastaout "drp_"${spl1}
done
for i in `seq 1 3`;
do
	cat drp_split_orf_* >drp${i}_all_orf.fa
	wc -l drp${i}_all_orf.fa
	rm drp_split_orf_*
	faToTab -type=protein drp${i}_all_orf.fa tmp
	sort -k 2 -d tmp |
		awk '{print ">"$1"\n"$2}'>drp${i}_all_orf.fa
	split -d -l 3000000 drp${i}_all_orf.fa split_orf_
		for spl1 in split_orf_*;
		do
			usearch10 -fastx_uniques ${spl1} -sizein -sizeout -minuniquesize 4 -threads 24 -fastaout "drp_"${spl1}
		done
	rm split_orf_* tmp*
done
cat drp_split_orf_* >${exp}"_drp.orf.fa"
rm all_drp_orf.fa
printf "%b\n" "-----------------------------------------------------------------------"
printf "%b\n" "Dereplicated ORF files generated:"
wc -l ${exp}"_drp.orf.fa"
printf "%b\n" "-----------------------------------------------------------------------"
printf '%(%a,%D;%T)T.\n' -1


mv ${exp}".clean_drp.fa" ${main_path}master/single_dereps/fullseqs/.
mv ${exp}"_drp.orf.fa" ${main_path}master/single_dereps/orf/.
rm drp_split* *_drp.orf *tmp*
mv *.fa results/.
mv *.orf results/.
wc -l results/*>LineCountResults.tsv
tar -zcvf TempReads.tar.gz TempReads/
rm -r TempReads
##Make sure that it goes back to the original folder at the end of each cycle
cd ${main_path}
