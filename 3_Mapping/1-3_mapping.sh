##  Pipeline for analysis of random peptide effect of e-coli in
##experiments described in Neme, 2017.
##Third step: Mapping orf extracted from the fastqfiles to DB
## (Script 3/)
##AUTHOR:JFajardo (01.2018)
##LAST_EDITED: 05.2020
##USAGE: ./1-3.orfMappingRP.sh [experiment folder number]

printf "%b\n" "-----------------------------------------------------------------------------"
printf "%b\n" "Mapping- Step 1: Individual mappings ${exp}"
printf "%b\n" "-----------------------------------------------------------------------------"
printf '%(%a,%D;%T)T.\n' -1

mkdir mapping
cd mapping/
awk '{print $1}' names_${exp}.tab > mapping/reps.tsv

## USEARCH (v.10) search to get count data for each experiment individually
#+using the -search_exact option (Available from v.8)
## Output format for the table is the same as blast6: 1=qseqid 2=sseqid 3=pident
#+ 4=length 5=mismatch 6=gapopen 7=qstart 8=qend 9=sstart 10=send 11=qlen
#+12=slen 11=evalue 12=bitscore

awk '{print $1}' BACT_DB_names.txt>names.tmp
cat names.tmp>table.tmp
while read -r rep;
do
	usearch10 -usearch_global ${rep}.merged.clean.fa -db BACT_DB.fa -threads 24 -id 0.97 -userout ${rep}.blast -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+ql+tl+evalue+bits+qcov+tcov -notmatched ${rep}.notmatched -strand plus -query_cov 0.9 -maxgaps 5
	printf "%b\n" "------------ Done: Mapping ${rep}.orf against database ----------------"

	awk -F "\t" '{print $2}' ${rep}.blast|
		cat - names.tmp|
		sort - >${rep}.tmp
	uniq -c ${rep}.tmp >${rep}.usearch.counts
	awk '{print $1}' ${rep}.usearch.counts >a_${rep}.tmp
	paste table.tmp a_${rep}.tmp > tmp2
	mv tmp2 table.tmp
	echo "end USEARCH: ${exp}, ${rep}"
	printf '%(%a,%D;%T)T.\n' -1
done < mapping/reps.tsv
mv table.tmp count_table.raw
rm *tmp tmp*

printf "%b\n" "Count of matched and notmatched queries"
wc -l *.blast *.notmatched

printf "%b\n" "-----------------------------------------------------------------------------"
printf "%b\n" "Mapping- Step 2: Count table generation"
printf "%b\n" "-----------------------------------------------------------------------------"
printf '%(%a,%D;%T)T.\n' -1
Rscript 1-3_mappingRP.R
mkdir countTables
mv count_table.tsv countTables/${exp}_count_table.tsv
echo "Done: Generation of USEARCH count table for "${exp}
printf '%(%a,%D;%T)T.\n' -1

## Run next script 1-4_ReportGen.sh
