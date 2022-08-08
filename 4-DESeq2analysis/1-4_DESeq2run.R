##AUTHOR: JFajardo (11.2020)
##LAST EDITED:11.2020
##CHANGES:Created

##Required packages: DESeq2, ggplot2, tidyverse
##DESeq2 analyses of amplicon sequencing data
## (Script 2) DESeq2 analysis

rm(list=ls(all=TRUE)) #Start with empty workspace

##Load necessary libraries
require("DESeq2")
require("ggplot2")
require("tidyverse")
require("Biostrings")
library("pheatmap")
library("vsn")
library("RColorBrewer")
library("scales")

#----Set working directories and source required functions----
home_loc<-c("")
exp_loc<-paste0(home_loc)
count_table_path<-paste0(exp_loc,"countTables/")

setwd(count_table_path)
getwd()
source("1-4_DESeq2functions.R")

##----Set output files and logs----
##Untreated samples
out_path<-paste0(count_table_path,"AnalysesDESeq2/")
check.make.dir(out_path)
#Create output folders beforehand
neg_folder<-check.make.dir(paste0(out_path,"Neg_tables/"))
pos_folder<-check.make.dir(paste0(out_path,"Pos_tables/"))
full_folder<-check.make.dir(paste0(out_path,"Full_tables/"))
sig_folder<-check.make.dir(paste0(out_path,"Sig_tables/"))
ns_folder<-check.make.dir(paste0(out_path,"NonSig_tables/"))

##Treated samples
#Create output folders beforehand
neg_folder2<-check.make.dir(paste0(out_path,"Neg_tables_tr/"))
pos_folder2<-check.make.dir(paste0(out_path,"Pos_tables_tr/"))
full_folder2<-check.make.dir(paste0(out_path,"Full_tables_tr/"))
sig_folder2<-check.make.dir(paste0(out_path,"Sig_tables_tr/"))
ns_folder2<-check.make.dir(paste0(out_path,"NonSig_tables_tr/"))

##----Read in databases----
##Complete database
my_db_full<-read_tsv("BACT_DB_w_names4.tsv")
my_db_full

##Prepare a summary table
summary_table<-read_tsv(file=paste0(exp_loc,"EXPinfo.tsv"))
summary_table<-tibble(summary_table,"min.5.count"=0, "min.5.count.per.db"=0,"signf.clones"=0,"signf.clones.per"=0,"up.clones"=0, "up.clones.per"=0,"down.clones"=0,"down.clones.per"=0,"ns.clones"=0,"ns.clones.per"=0)
summary_table

##Get list of files in the count table folder
my_tables<-list.files(pattern = "count_table.tsv")
files.to.use<-pmatch(paste0(summary_table$folder.name,"_"), my_tables)
my_tables<-my_tables[files.to.use]
my_tables

# ##Use this for testing
table_name<-my_tables[8]
table_name

##Prepare a variable to store samples with treatments
treated_tables<-vector(mode = "character")
treated_tables

##Run DESeq2 in non-treated samples
for(table_name in my_tables){
  exp_prefix0<-unlist(str_split(table_name,pattern = "_"))[1]
  exp_prefix<-gsub("-","",exp_prefix0)
  exp_prefix

  ##Set Output file
  sink(file=paste0(out_path,"DESeq2log",exp_prefix,".txt"))
  # sink()
  ##Output figures
  pdf(file = paste0(out_path,"2020_",exp_prefix,"DESeq2-%03d.pdf"), family = "Helvetica", height = 5, width = 7, useDingbats = FALSE, onefile = FALSE)
  # svg(filename = paste0(out_path,"2020_",exp_prefix,"DESeq2%03d.svg"), family = "Helvetica", height=5, width = 7)

  ##----Prepare data----
  ##Read in count data (in count table format), retrieve coldata from column names
  ##(COLUMNS HAVE the following naming format:
  ##Cycle-Replicate-Treatment_orf_count_table.dat), prepare dataset for analysis.
  countdata<-read.table(table_name,row.names = 1,stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

  ##Remove sequences that have less than 5 counts
  countdata <- countdata[rowSums(countdata)>5,]

  ##Remove empty plasmid (Outlier)
 # countdata <- countdata[-1,]

  ##Remove plasmid sample (no replicates)
  if("plasmid"%in%colnames(countdata)==1){
    countdata <- countdata[,-which(colnames(countdata)=="plasmid")]
  }else if("0-1"%in%colnames(countdata)==1){
    countdata <- countdata[,-which(colnames(countdata)=="0-1")]
  }else if("8am"%in%colnames(countdata)==1){
    countdata <- countdata[,-which(colnames(countdata)=="8am")]
  }
  countdata
  dim(countdata)

  ##Add experiment name to column names and extract column info from the column names
  colnames(countdata)<-paste0(exp_prefix,"-",colnames(countdata))
  head(countdata)
  replist<-colnames(countdata)
  replist
  cdat0<-strsplit(replist,split="-")
  cdat0
  coldata <- cbind(unlist(replist), unlist(lapply(cdat0, function(x) x[1])), unlist(lapply(cdat0, function(x) x[2])), unlist(lapply(cdat0, function(x) x[3])), unlist(lapply(cdat0, function(x) x[4])))
  coldata
  colnames(coldata)<-c("Sample", "Experiment", "Cycle", "Rep", "Treatment")
  coldata<-as.data.frame(coldata)
  coldata$Treatment<-as.factor(coldata$Treatment)
  coldata

  ##Modify column names in count table to match information in coldata table
  colnames(countdata)<-coldata$Sample
  rownames(coldata)<-coldata$Sample
  head(countdata)
  str(countdata)

  ##Save table name in treated variable if present
  treatments<-levels(coldata$Treatment)
  treatments
  if(length(treatments)>0){
    treated_tables<-c(treated_tables,table_name)
    paste("Table added to treatment table list:",treated_tables)
    }else{
    print("Untreated samples only")
  }

  ##Generate DESeq data set from countdata
  dds<- DESeqDataSetFromMatrix(
    countData = countdata,
    colData = coldata,
    design = ~Cycle)
  dds

  ##Make sure that "Control" is the first level in the Cycle
  if ("8am"%in%dds$Cycle){
    dds$Cycle<- relevel( dds$Cycle, "11am")
    print("Releveled 11am to first level")
  } else if ("1"%in%dds$Cycle){
    dds$Cycle<- relevel( dds$Cycle, "1")
    print("Releveled 1 to first level")
  } else if ("plasmid"%in%dds$Cycle){
    dds$Cycle<- relevel( dds$Cycle, "plasmid")
    print("Releveled plasmid to first level")
  } else if ("0-1"%in%dds$Cycle){
    dds$Cycle<- relevel( dds$Cycle, "0-1")
    print("Releveled 0-1 to first level")
  } else {
    print("Levels not modified for dds")
  }

  ##Subset relevant columns from the data matrix. In this case, only untreated samples
  dds <- dds[,is.na(dds$Treatment)]
  dds <- dds[,!is.na(dds$Rep)]
  #Remove factor levels of the removed columns
  dds$Cycle <- droplevels(dds$Cycle)
  dds
  print(dds$Cycle)
  ##Get cycle levels
  cycles <- levels(dds$Cycle)
  cycles

  ##Run DESeq2 for the untreated subset of samples
  print("Run DESeq2 for the untreated subset of samples")
  dds <- DESeq(dds)
  res.names<-resultsNames(dds)
  coef.name<-res.names[length(res.names)]
  res <- results(dds, name = coef.name)
  summary(res)
  resLFC <- lfcShrink(dds, coef = coef.name, type="apeglm")
  plotMA(resLFC, ylim=c(-2,2))
  plotMA(res, ylim=c(-2,2))
  plotCounts(dds, gene=which.min(res$padj), intgroup="Cycle")
  plotCounts(dds, gene = "BACT000000001", intgroup = "Cycle")
  vsd <- vst(dds, blind=FALSE)
  # rld <- rlog(dds, blind=FALSE)
  ntd <- normTransform(dds)
  meanSdPlot(assay(ntd))
  meanSdPlot(assay(vsd))
  # meanSdPlot(assay(rld))
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[,c("Cycle", "Rep")])
  pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
           cluster_cols=TRUE, annotation_col=df)
  pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
           cluster_cols=TRUE, annotation_col=df)

  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$Cycle, vsd$Rep, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  plotPCA(vsd, intgroup=c("Cycle", "Rep"))

  ##Set up individual contrasts between different Cycles
  ##set.contratst function saved in 1-6_DESeq2functions.R
  print("Set-up contrasts")
  contr <- set.contrasts(dds)

    ##Get Contrast names to use as plot labels
  contr_info <- lapply(contr, function(x) mcols(x, use.names = TRUE))
  contr_info
  my_legends <- sapply(contr_info, function(x) unlist(strsplit(x[2,2], ": "))[2])
  my_legends
  my_prefix<-gsub(" ","",my_legends)
  my_prefix
  names(contr) <- my_legends

  if ("Intercept" %in% names(contr)){
    print("removing Intercept from contrast list")
    contr<-contr[-which(names(contr)=="Intercept")]
  }
  names(contr)

  ##Add a column describing the magnitude of adjusted pvalue, and direction of fold change to the contrast table. (type.padj function in function file)
  contr<-lapply(contr,function(x) data.frame(x,"colVector" = type.padj(x)))
  contr<-lapply(contr,function(x) data.frame(x,"colVector2" = type.padj(my.matrix = x, padj = 1, non.sig = "NS", pos.code = "UP", neg.code = "DOWN")))
  tail(contr)
  contr<-lapply(contr, function (x) data.frame(x,my_db_full[match(rownames(x), my_db_full$seq.name),3:16]))
  contr
  contr<-lapply(contr, function(x) as_tibble(rownames_to_column(x,var = "seq.name")))
  contr

  ##Save tables of different significance levels
  contr05 <- lapply(contr,function(x) filter.contr(x,padj<0.05))
  contr05
  contr01 <- lapply(contr,function(x) filter.contr(x,padj<0.01))
  contr01
  contrPos05 <- lapply(contr05,function(x) filter.contr(x,log2FoldChange> 1))
  contrPos05
  contrNeg05 <- lapply(contr05,function(x) filter.contr(x,log2FoldChange< (-1)))
  contrNeg05
  contrPos05HFC <- lapply(contr05,function(x) filter.contr(x,log2FoldChange>2))
  contrPos05HFC
  contrNeg05HFC <- lapply(contr05,function(x) filter.contr(x,log2FoldChange<(-2)))
  contrNeg05HFC
  contrPos01 <- lapply(contr01,function(x) filter.contr(x,log2FoldChange> 1))
  contrPos01
  contrNeg01 <- lapply(contr01,function(x) filter.contr(x,log2FoldChange< (-1)))
  contrNeg01
  contrPos01HFC <- lapply(contr01,function(x) filter.contr(x,log2FoldChange>2))
  contrPos01HFC
  contrNeg01HFC <- lapply(contr01,function(x) filter.contr(x,log2FoldChange<(-2)))
  contrNeg01HFC
  contrNS <- lapply(contr,function(x) filter.contr(x,padj>0.05))
  contrNS
  contrNSlowSD <- lapply(contr,function(x) filter.contr(x,log2FoldChange> -1&log2FoldChange<1&padj>0.05))
  contrNSlowSD

  ##Generate color-coded MA plots
  for (m in 1:length(contr)){
    my.MAplot(contr[[m]], my.title=paste(exp_prefix,names(contr[m])) ,peval=0.05)
    my.MAplot(contr[[m]], my.title=paste(exp_prefix,names(contr[m])) ,peval=0.01)
    my.MAplotLen(contr[[m]], my.title=paste(exp_prefix,names(contr[m])) ,peval=0.05)
  }
  ##Plot individual clones
  head(res[order(res$padj),], 10)
  sum(res$padj<0.05, na.rm = TRUE)
  toplot<-rownames(head(res[order(res$padj),], 25))
  topGenes <- head(order(res$padj),25)
  diffFull<-rownames(res[which(res$padj<0.05),])
  # Plot trajectories of top 25 genes
  cnts <- counts(dds, normalized = TRUE, replaced = FALSE)[toplot, ]
  cnts<-as_tibble(t(cnts))
  cnts
  Cycle<-colData(dds)[["Cycle"]]
  cnts<-bind_cols(cnts, Cycle=Cycle)%>%
    pivot_longer(cols = 1:25, names_to = "Clone", values_to = "Counts")
  cnts$Cycle<-as.factor(cnts$Cycle)
  p<-ggplot(cnts, aes(x = Cycle, y = Counts, group=1)) +
    geom_point()  +
    stat_summary(fun=mean, geom="line", color="blue")+
    theme_minimal(base_size = 10, base_family = "Helvetica")+
    facet_wrap(facets = vars(Clone), scales = "free_y")
  print(p)
  p
  ##Name the columns of the contrast lists to their values and write the lists into individual tables
  write.list2tables(write.path=full_folder, contr, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq.tsv")
  write.list2tables(write.path=sig_folder, contr05, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_Sig05.tsv")
  write.list2tables(write.path=pos_folder, contrPos05, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_Pos05.tsv")
  write.list2tables(write.path=neg_folder, contrNeg05, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_Neg05.tsv")
  write.list2tables(write.path=pos_folder, contrPos05HFC, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_Pos05HFC.tsv")
  write.list2tables(write.path=neg_folder, contrNeg05HFC, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_Neg05HFC.tsv")
  write.list2tables(write.path=sig_folder, contr01, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_Sig01.tsv")
  write.list2tables(write.path=pos_folder, contrPos01, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_Pos01.tsv")
  write.list2tables(write.path=neg_folder, contrNeg01, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_Neg01.tsv")
  write.list2tables(write.path=pos_folder, contrPos01HFC, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_Pos01HFC.tsv")
  write.list2tables(write.path=neg_folder, contrNeg01HFC, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_Neg01HFC.tsv")
  write.list2tables(write.path=ns_folder, contrNS, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_NonSig.tsv")
  write.list2tables(write.path=ns_folder, contrNSlowSD, tab.name = my_prefix, file.prefix = exp_prefix,file.ext = "_DESeq_NonSig_lowSD.tsv")

  ##Add information the the summary table
  rowindex<-which(summary_table$folder.name==exp_prefix0)

  summary_table[rowindex,c(9:18)]<- list(nrow(res), nrow(res)*100/nrow(my_db_full), nrow(contr05[[my_legends[length(my_legends)]]]),  nrow(contr05[[my_legends[length(my_legends)]]])*100/nrow(res), nrow(contrPos05[[my_legends[length(my_legends)]]]), nrow(contrPos05[[my_legends[length(my_legends)]]])*100/nrow(res),  nrow(contrNeg05[[my_legends[length(my_legends)]]]), nrow(contrNeg05[[my_legends[length(my_legends)]]])*100/nrow(res), nrow(contrNS[[my_legends[length(my_legends)]]]), nrow(contrNS[[my_legends[length(my_legends)]]])*100/nrow(res))

  ##Add colvector to generaldb table
  colVectib<-contr[[my_legends[length(my_legends)]]][c("seq.name","colVector","colVector2")]
  names(colVectib)<-c("seq.name",paste0(exp_prefix,".0.05sign"),paste0(exp_prefix,".trend"))
  my_db_full<-my_db_full%>%
    left_join(.,colVectib, by="seq.name")


  graphics.off()
  sink()
}
graphics.off()
sink()
write_tsv(summary_table, file = paste0(out_path,"ALL_exp_summary.tsv"))
write_tsv(my_db_full, file = "BACT_DB_w_names5.tsv")

##Run next script: 1-6_DB_extras.sh
