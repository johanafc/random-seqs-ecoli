##  Pipeline for analysis of random peptide effect of e-coli in
#experiments described in Neme, 2017.
##Generating a database of peptides
## Generating count data table of peptides from new gene data of e-coli in
#experiments described in Neme, 2017.

#Input data is a table with no column names and the count data for each ref sequence.
#Output is the count table with the Sequence names as row names
#and the column names with the conditions

#Read data table
my_table<-read.table(file = "count_table.raw",header = FALSE)

# The table generation artificially adds 1 to every sequence
my_table[,2:ncol(my_table)]<-my_table[,2:ncol(my_table)]-1

#Set colnames as day and replicate (d1-r1) for each experiment
my_names<-read.table("reps.tsv", stringsAsFactors = FALSE)
my_names<-c("SeqName",my_names[[1]])
colnames(my_table) <- my_names
write.table(my_table, file="count_table.tsv", sep="\t", quote=F, row.names=FALSE)
