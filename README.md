# Pipeline used for the analysis of amplicon sequencing data generated from a library of random sequences in _E. coli_

The scripts in this repository were used for the analysis of amplicon sequencing date as described in "Castro JF, Tautz D. The Effects of Sequence Length and Composition of Random Sequence Peptides on the Growth of E. coli Cells. Genes (Basel). 2021 Nov 28;12(12):1913. doi: 10.3390/genes12121913. PMID: 34946861; PMCID: PMC8702183."

Available data for separate experiments was saved in individual folders labelled "exp-X" with "X" being an integer form 1 to 9. 

The pipeline consists of the following steps:

1. Extraction of random sequences between the sequencing primers.
2. Generation of the database of unique sequences present in the sequencing data.
3. Description of sequence features (length, GC content, ISD, for each sequence in the database).
4. Mapping the sequencing data to the database to get count tables.
5. Calculating frequency changes for each sequence in the database.
6. Cross-comparison of the results obtained for each individual experiment.

Dependencies:
This pipeline uses several third party software that should be installed and be functional:
- Trimmomatic (v 0.36)
- USEARCH10
- getorf
- transeq
- faToTab
- seqtk
- IUpred2

Additionallly, protein aggregation propensities in step 2 were calculated using an online tool: PASTA2.0 http://old.protein.bio.unipd.it/pasta2/ (Ian Walsh, Flavio Seno, Silvio C.E. Tosatto and Antonio Trovato. PASTA2: An improved server for protein aggregation prediction. Nucleic Acids Research, 2014 Jul;42(Web Server issue):W301-7.)
