require("ggplot2")
require("tidyverse")
require("Biostrings")
library("pheatmap")
library("vsn")
library("RColorBrewer")
library("scales")

check.outdir <- function(my.path=getwd(),def.out) {
  if (def.out!=""){
    if(!file.exists(def.out)){
      cat("Creating user defined Output folder: ", def.out)
      out.path<-def.out
      dir.create(out.path, recursive=TRUE)
    }else{
      cat("User defined Output folder:", def.out)
      out.path<-def.out
    }
  }else if (!file.exists(paste0(my.path,"Results/"))){
    out.path<-paste0(my.path,"Results/")
    cat("Creating output folder", out.path)
    dir.create(out.path, recursive=TRUE)
  }else{
    out.path<-paste0(my.path,"Results/")
    cat("Using existing output folder:", out.path)
  }
  print(out.path)
}

check.make.dir <- function(def.folder) {
  if(!file.exists(def.folder)){
    cat("Creating user defined folder: ", def.folder)
    dir.create(def.folder, recursive = TRUE)
  }else{
    cat("User defined folder exists:", def.folder)
  }
  print(def.folder)
}

my.read.table <- function(my.count.table,...) {
  read.table(my.count.table,row.names=1,check.names=F,stringsAsFactors=F,sep = "\t",...)
}

get.coldata <- function(hysep.names,str.sep="-",i=c(1:4)) {
  cdat0<-strsplit(hysep.names,split=str.sep)
  cdat <- cbind(unlist(hysep.names), unlist(lapply(cdat0, function(x) x[1])), unlist(lapply(cdat0, function(x) x[2])), unlist(lapply(cdat0, function(x) x[3])), unlist(lapply(cdat0, function(x) x[4])))
  print(cdat)
}

set.contrasts <- function(dds.table) {
  require("DESeq2")
  cont.names <- resultsNames(dds.table)
  cont.res<-vector(mode = "list",length = length(cont.names)) 
  i<-1
  for (l in cont.names){
    cont.res[[i]]<- results(dds.table, name = l)
    cont.res[[i]] <- cont.res[[i]][order(cont.res[[i]]$log2FoldChange),]
    i<-i+1
  }
  return(cont.res)
}

type.padj<-function(my.matrix,padj=0.05,non.sig="ns",pos.code="pos",neg.code="neg"){
  colVector<-vector(mode="character",length = length(my.matrix[,"padj"]))
  #  print(colVector)
  colVector[which(is.na(my.matrix[,"padj"]))]<-non.sig
  colVector[which(my.matrix[,"padj"]<padj&my.matrix[,"log2FoldChange"]>1)]<-pos.code
  colVector[which(my.matrix[,"padj"]<padj&my.matrix[,"log2FoldChange"]< (-1))]<-neg.code
  colVector[which(colVector=="")]<-non.sig
  return(factor(colVector,levels = c(pos.code,neg.code,non.sig)))
}

filter.contr<-function(contr.tibble,filter.cond){
  filter.cond<-enquo(filter.cond)
  contr.tibble%>%
    filter(.,!!filter.cond)
}

write.list2tables <- function(write.path=out.path1,cont.tab,tab.name="",file.prefix="",file.ext="tsv") {
  for (x in 1:length(cont.tab)) {
    tab.pref<-tab.name[x]
    write_tsv(cont.tab[[x]], file = paste0(write.path,file.prefix,tab.pref,file.ext), col_names = TRUE)
  }
}

my.MAplot <- function(deseq.results.table, my.title="", my.colours =c("mediumorchid4","deepskyblue3", "gray27"),peval=0.05,...) {
  gplot<-ggplot(deseq.results.table, aes(x=log2FoldChange, y=baseMean))+
    geom_point(aes(col=colVector, shape=colVector),size=1)+
    labs(title=my.title, y="Base mean counts", x="Log2 fold-change", caption=paste(length(which(deseq.results.table$padj<peval)),"sign. \n", length(which(deseq.results.table$padj<peval&deseq.results.table$log2FoldChange>1)),"pos. \n", length(which(deseq.results.table$padj<peval&deseq.results.table$log2FoldChange< (-1))),"neg. \n",length(deseq.results.table$padj),paste0("total \n padj<",peval)), color="Group", shape="Group")+
    scale_colour_manual(values = my.colours, labels=c("pos"="Positive, significant", "neg"="Negative, significant", "ns"="Non-Significant"))+
    scale_shape_manual(values = c("pos"=2, "neg"=6, "ns"=20),labels=c( "pos"="Positive, significant", "neg"="Negative, significant", "ns"="Non-Significant"))+
    scale_y_log10(labels = scientific)+
    coord_cartesian(xlim = c(-8,8))+
    theme_linedraw()
  print(gplot)
}

my.MAplotLen <- function(deseq.results.table, my.title="", low.colour="lightgreen", high.colour="chocolate2", peval=0.05,...) {
  gplot2<-ggplot(deseq.results.table, aes(x=log2FoldChange, y=baseMean))+
    geom_jitter(aes(col=pep.len, shape=colVector),size=1.2, alpha=0.8)+
    labs(title=my.title, y="Log10 base mean counts", x="Log2 fold-change", caption=paste(length(which(deseq.results.table$padj<peval)),"sign. \n", length(which(deseq.results.table$padj<peval&deseq.results.table$log2FoldChange>1)),"pos. \n", length(which(deseq.results.table$padj<peval&deseq.results.table$log2FoldChange< (-1))),"neg. \n",length(deseq.results.table$padj),paste0("total \n padj<",peval)), color="Group", shape="Group")+
    scale_color_gradient(low = low.colour, high = high.colour)+
    scale_shape_manual(values = c("pos"=2, "neg"=6, "ns"=20),labels=c( "pos"="Positive, significant", "neg"="Negative, significant", "ns"="Non-Significant"))+
    scale_y_log10(labels = scientific)+
    coord_cartesian(xlim = c(-8,8))+
    theme_linedraw()
  print(gplot2)
}

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


plotLen<-function(my.data, my.title="", my.ylim=c(0,800)){
  geofun<-function(x){
    nrow(my.data)*(1-3/64)^(x-1)*3/64
  }
  expfun<-function(x){
    (20^(x+1))
  }
  ggplot(data = my.data)+
    xlim(0,67)+
    geom_histogram(aes(x=pep.len), binwidth=1, colour="black", fill="lightblue")+
    stat_function(fun=geofun,geom="line", color="darkgrey", linetype="solid", size=0.5)+
    stat_function(fun=expfun, color="black", linetype="dashed")+
    labs(title = my.title, x="Length (AA)", y="Count",caption = paste0("Total number of sequences=",nrow(my.data)))+
    coord_cartesian(ylim = my.ylim)+
    theme_minimal(base_size = 13, base_family = "Helvetica")
}


plotGC<-function(my.data,my.title=""){
  ggplot(data = my.data)+
    geom_histogram(aes(x=per.gc.full, fill="Full_Reads"), binwidth=1, alpha=0.2, colour="black")+
    geom_histogram(aes(x=per.gc.orf, fill="ORFs"), binwidth=1,alpha=0.2,colour="black")+
    geom_vline(xintercept=50,linetype="dashed", colour="seagreen")+
    labs(title = my.title, x="GC%", y="Count", fill="")+
    coord_cartesian(xlim=c(20,80))+
    theme_light(base_size = 13, base_family = "Helvetica")
}

plotIDS<-function(my.data, my.title = ""){
  ggplot(data = my.data)+
    geom_histogram(aes(x=iupred.long.frac), binwidth=0.01, colour="black", fill="grey")+
    labs(title = my.title, x="Average IDS-short", y="Count",caption = paste0("Total number of sequences=",nrow(my.data)))+
    coord_cartesian()+
    theme_minimal(base_size = 13, base_family = "Helvetica")
}
