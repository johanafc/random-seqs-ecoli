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
