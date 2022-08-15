#!bin/bash
##  Pipeline for analysis of random peptide effect of e-coli in
#+ experiments described in Neme, 2017.
##Second step: Generate a master DB with data from all combined experiments
##AUTHOR: Johana Fajardo (10.2017)
##LAST EDIT: 17.07.2020
##Output files:

main_path="" #Path to the folder where all experiments are located
db_path="Databases/" #Path to the folder where databases will be located


## Get a dereplicated file of all orf sequences from all experiments
## USEARCH10: 32bit version of USEARCH10 only allows to use up to 4Gb of memory
#+ at the time, if one tries to dereplicate a bigger file, it fails. To prevent
#+ this, we have to split the files in 3K lines, dereplicate each split file and
#+ merge.
printf "%b\n" "-----------------------------------------------------------------------------"
printf "%b\n" "DB generation- Step 1: Dereplicating"
printf "%b\n" "-----------------------------------------------------------------------------"
printf '%(%a,%D;%T)T.\n' -1

cat *".clean_drp.fa" > all_drp.fa
faToTab -type=protein all_drp.fa tmp1
sort -k 2 -d tmp1 > tmp1.tsv
awk '{print ">"$1"\n"$2}' tmp1.tsv >all_drp_sorted.fa
split -d -l 3000000 all_drp_sorted.fa split_drp_
##9 splits (min 1 median 1 max>59000 avg>10)
for spl in split_drp_*;
do
	usearch10 -fastx_uniques ${spl} -sizein -sizeout -threads 24 -minuniquesize 8 -fastaout "drp_"${spl}
done
printf "%b\n" "Count lines splits and de-replicated splits"
wc -l split_drp_* drp_split_drp_*
rm split_drp_* tmp*

for i in `seq 1 3`;
do
	cat drp_split_drp_* > drp${i}_all_drp.fa
	rm drp_split_drp_*
	faToTab -type=protein drp${i}_all_drp.fa tmp2
	sort -k 2 -d tmp2>tmp2.sort
	awk '{print ">"$1"\n"$2}' tmp2.sort >tmp2.fa
	split -d -l 3000000 tmp2.fa split_drp_
	for spl1 in `ls -1v split_drp_*`;
	do
		usearch10 -fastx_uniques ${spl1} -sizein -sizeout -threads 24 -minuniquesize 8 -fastaout "drp_"${spl1}
	done
	printf "%b\n" "Count lines splits and de-replicated splits (${i})"
	wc -l split_drp_* drp_split_drp_*
	rm split_drp_* tmp*
done
cat drp_split_drp_* > tmp
seqtk seq tmp>full_exact_drp.fa

faToTab -type=protein full_exact_drp.fa full_exact_drp.tsv
printf "%b\n" "Count lines final exact match de-replicated DB"
wc -l full_exact_drp.tsv
awk '{print length($2)}' full_exact_drp.tsv> full_exact_drp_lens.tsv
rm drp_split_drp_* tmp

usearch -sortbysize full_exact_drp.fa -fastaout full_exact_drp_sorted.fa
usearch -sortbysize full_exact_drp.fa -fastaout full_exact_drp_sorted_noempty.fa -minsize 8 -maxsize 1000000

for myfile in full_exact_drp_sorted*.fa;
do
	i2=${myfile##*sorted}
	i1=${i2%%.fa}
	usearch -cluster_smallmem ${myfile} -id 0.97 -sizein -sizeout -sortedby size -centroids tmp -uc full_97ID_clusters${i1}.uc -minsize 15
	seqtk seq tmp >full_97ID_clustered${i1}.fa
done
wc full_97ID_clustered*

for myfile in full_exact_drp_sorted*.fa;
do
	i2=${myfile##*sorted}
	i1=${i2%%.fa}
	usearch -cluster_otus ${myfile} -minsize 15 -sizein -sizeout -otus tmp -uparseout full_uparse${i1}.txt
	seqtk seq tmp >full_97ID_otus${i1}.fa
done
wc full_97ID_otus*

for myfile in full_97ID_*.fa; do faToTab -type=protein $myfile tmp; mv tmp ${myfile%.*}.tsv ; done
for myfile in full_97ID_*.tsv; do awk '{print length($2)}' $myfile >tmp; mv tmp ${myfile%.*}_lens.tsv ; done
printf "%b\n" "Finished second round of dereplication"

## Add the prefix "BACTORFNR" to every sequence name, save in BACT_ORF_NR_DB.fa
printf "%b\n" "-----------------------------------------------------------------------------"
printf "%b\n" "DB generation- Step 2: Generating DB fasta, sequence, and name files"
printf "%b\n" "-----------------------------------------------------------------------------"
printf '%(%a,%D;%T)T.\n' -1
awk 'BEGIN{a=1000000000}{a=a+1; print "BACT"a"\t"$1"\t"$2}' full_97ID_clustered.tsv |
	sed 's/BACT1/BACT/' >BACT_DB_w_names0.tsv
printf "%b\n" "FULL DB in tabular format saved as: BACT_DB_w_names0.tsv"
awk '{print ">"$1"\n"$3}' BACT_DB_w_names0.tsv >BACT_DB.fa
printf "%b\n" "FULL DB in fasta format saved as: BACT_DB.fa"
awk '{print $3}' BACT_DB_w_names0.tsv>BACT_DB_seqsOnly.txt
printf "%b\n" "FULL DB sequence list saved as: BACT_DB_seqsOnly.txt"
awk '{print $1"\t"$2}' BACT_DB_w_names0.tsv>BACT_DB_names.txt
printf "%b\n" "FULL DB name list saved as: BACT_DB_names.txt"
awk '{print length($3)}' BACT_DB_w_names0.tsv>BACT_DB_readLens.txt
printf "%b\n" "Read length list saved as: BACT_DB_readLens.txt"
awk '{gsub("B.*size=",""); gsub(";.*",""); print}' BACT_DB_w_names0.tsv > BACT_DB_clusterSizes.txt
printf "%b\n" "Cluster size list saved as: BACT_DB_clusterSizes.txt"
cusp -sequence BACT_DB.fa -outfile BACT_DB.cut
printf "%b\n" "Codon usage table of full reads saved as: BACT_DB.cut"


printf "%b\n" "-----------------------------------------------------------------------"
printf "%b\n" "Getting random peptide ORFs from Dereplicated database"
printf "%b\n" "-----------------------------------------------------------------------"
printf '%(%a,%D;%T)T.\n' -1

#Generate fasta and tsv for all ORFs present in the database
getorf -sequence BACT_DB.fa -outseq tmp -minsize 12 -find 3
seqtk seq tmp> BACT_DB_allORFs.fa
faToTab -type=protein BACT_DB_allORFs.fa BACT_DB_allORFs.tsv

#Find ORFs that correspond to the expected sequences in the DB. Generate an error if any sequence in the DB does not have them.

rm -v BACT_ORF_DB.tsv tmp
touch BACT_ORF_DB.tsv

while read -r line;
do
	seqname=${line%%_*}
	if ! grep -q $seqname BACT_ORF_DB.tsv;
	start='[[:space:]]ATGAAGCTT';
	then
		if [[  "$line" =~ $start ]];
		then
			echo $line>>BACT_ORF_DB.tsv;
		elif [[ "$line" == *"_1"* ]];
		then
			echo $line>>BACT_ORF_DB.tsv;
		else
			echo "Unexpected error. No valid ORF found"
#			echo "$seqname\tNNN">>BACT_ORF_DB.tsv
		fi
	fi
done <BACT_DB_allORFs.tsv
echo "Last sequence in file: $seqname"

rm -v *tmp*

wc -l BACT_ORF_DB.tsv
wc -l BACT_DB_w_names0.tsv
echo "Check here whether the line numbers are the same\n\n"

#Get fasta file of the clean ORF DB, and a list of ORF sequences only
awk '{print ">"$1"\n"$2}' BACT_ORF_DB.tsv>BACT_ORF_DB.fa
cusp -sequence BACT_ORF_DB.fa -outfile BACT_ORF_DB.cut

paste BACT_DB_names.txt BACT_DB_clusterSizes.txt BACT_DB_seqsOnly.txt BACT_DB_readLens.txt >BACT_DB_w_names1.tsv

awk '{print $2}' BACT_ORF_DB.tsv>tmp
paste BACT_DB_w_names1.tsv tmp> tmp2
#File BACT_DB_w_names2.tsv includes ORFs and lengths
awk '{print $0"\t"length($6)}' tmp2>BACT_DB_w_names2.tsv

#Get predicted ORFs sorted alfabetically by sequence
sort -k 6 -d BACT_DB_w_names2.tsv|
 awk '{print ">"$1"\n"$6}'>predictedORF_seqsorted.fa

sort -k 1 -d BACT_DB_w_names2.tsv|
	awk '{print $6}'>BACT_ORF_DB_seqsOnly.txt

rm -v *tmp*

#Generate de-replicated ORF database (removes reads that code for the same peptide sequence, keeping track of how many)
usearch10 -fastx_uniques predictedORF_seqsorted.fa -sizeout -threads 24 -fastaout tmp -tabbedout BACT_ORF_DB_drp_uc.tsv

faToTab -type=protein tmp tmp2
sort tmp2|
	sed 's/BACT/BACTORF/'|
	awk '{print ">"$1"\n"$2}'>BACT_ORF_DB_drp.fa
cusp -sequence BACT_ORF_DB_drp.fa -outfile BACT_ORF_DB_drp.cut


faToTab -type=protein BACT_ORF_DB_drp.fa BACT_ORF_DB_drp.tsv
wc -l *
rm -v *tmp*
printf "%b\n" "Dereplicated ORF database generated"

## Translate the orf database and dereplicate it. Change the prefix of the seqs
#+ for "PEPNR" and save it as the database BACT_PEP_DB_drp.fa
printf "%b\n" "-----------------------------------------------------------------------------"
printf "%b\n" "Generating Peptide DB fasta, sequence, and name files"
printf "%b\n" "-----------------------------------------------------------------------------"
printf '%(%a,%D;%T)T.\n' -1

transeq -sequence BACT_ORF_DB.fa -outseq stdout -trim TRUE -supper1 -frame 1 -table 11|
	seqtk seq >BACT_PEP_DB.fa

sed '/>/d' BACT_PEP_DB.fa>BACT_PEP_DB_seqsOnly.txt

paste BACT_DB_w_names2.tsv BACT_PEP_DB_seqsOnly.txt> tmp
awk '{print $0"\t"length($8)}' tmp>BACT_DB_w_names3.tsv

sort -k 6 -d BACT_DB_w_names3.tsv|
 awk '{print ">"$1"\n"$8}'>predictedPEP_seqsorted.fa

#Save Dereplicated peptide database to its own fasta file
usearch10 -fastx_uniques predictedPEP_seqsorted.fa -sizein -sizeout -threads 24 -fastaout tmp3 -tabbedout BACT_PEP_DB_drp_uc.tsv

faToTab -type=protein tmp3 tmp23
sort tmp23|
	sed 's/BACT/BACTPEP/'|
	awk '{print ">"$1"\n"$2}'>BACT_PEP_DB_drp.fa

printf "%b\n" "Peptide DB in fasta format saved as: BACT_PEP_DB_drp.fa"

## Add cluster sizes to the table



mv *_DB* ${db_path}.
mv predicted* ${db_path}.
rm *tmp*

#Get iupred short values from peptide sequences
##Check that column name indicated matches the peptide
#sequence column on the table
bash 1-2_iupred.sh BACT_DB_w_names3.tsv

#Get GC content full sequences
perl 1-2_GCcontent.pl BACT_DB_seqsOnly.txt full
#Get GC content ORFs
perl 1-2_GCcontent.pl BACT_ORF_DB_seqsOnly.txt orf

##Get aggregation scores from PASTA2.0 server (upload PEP
##fasta)
cd PASTA2
##First make sure that all results have only one best pairing energy
for myfile in *.seq.best_pairings_list.dat; do wc -l $myfile| grep -v ^" *1" -; done
## No output. Means there is only one line in each file
## Best pairing files look like this:
## pairing 0  PASTA energy -2.457691  length 6  between segments 3-8 and 8-3  antiparallel

awk '{print $1}' BACT_DB_names.txt>DB_names_only.txt
rm best_pairings.txt
while read -r line;
do
    awk '{print $5}' BACT_PEP_DB_PASTApred/${line}.fasta.seq.best_pairings_list.dat >>best_pairings.txt
done <BACT_DB_names_only.txt

#Ian Walsh, Flavio Seno, Silvio C.E. Tosatto and Antonio Trovato.PASTA2: An improved server for protein aggregation prediction. Nucleic Acids Research, accepted. (2014)

paste BACT_DB_w_names3.tsv ids-short-avg ids-short-frac ids-long-avg ids-long-frac GCcontent-full.fa GCcontent-orf.fa best_pairings.txt>BACT_DB_w_names4.tsv

rm -v GCcontent*.fa ids-* best_pairings.txt

## Add column names to the complete table
##

sed -i -e '1 i seq.name\tread.name\tcluster.size\tfull.seq\tfull.len\torf.seq\torf.len\tpep.seq\tpep.len\tiupred.short.avg\tiupred.short.frac\tiupred.long.avg\tiupred.long.frac\tper.gc.full\tper.gc.orf\tpasta.best.e' BACT_DB_w_names4.tsv

##Run next script 1-3.mapping.sh
