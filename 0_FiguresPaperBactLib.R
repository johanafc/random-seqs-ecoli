## Created by: Johana Fajardo (07.2021)
## Last modified: (17.08.2021)
## Changes: Changed length distribution and aa freq plots.
##          Replaced Mann Whithney U test with KS test
## R version 4.0.3
##Un-comment lines starting with pdf and dev.off to save PDF files

## Start with a clean workspace
rm(list=ls(all=TRUE))

## Load required libraries
require(scales)
require(Biostrings)
require(pheatmap)
require(RColorBrewer)
require(ggpubr)
require(cowplot)
require(gridExtra)
require(dplyr)
require(tidyverse)

## Set working directory
setwd("/Users/jfajardo/ownCloud/2021-Paper-BL/")
# setwd("/home/jfajardo/ownCloud/2021-Paper-BL/")
getwd()

##Set color scales
colscale<-c("ORFs"="#89112d","POS"= "#f7a707","Simulation"= "#f5d36c","Low"="#f5e8c1","NEG"="#3ecad7","FLAG"="#487de0","RandomSeq"="#1b69a4", "Database"="#57359c","Black"="#000000", "NS"="#bcbdc1")
colscale0<-c("#f7a707","#f5d36c","#f5e8c1","#3ecad7","#44B6D0","#2d52d8","#57359c","#6F145D","#000000","#bcbdc1")
colscalelen<-colscale0[c(5:9)]
colscalegroup<-colscale0[c(1,10,4)]
colscaleflag<-colscale0[c(5,9)]

##Order of amino acids according to Uversky,2010 (order promoting to disorder promoting)
myorder<-c("W", "F", "Y", "I", "M", "L", "V", "N", "C", "T", "A", "G", "R", "D", "H", "Q", "K", "S", "E", "P")

##Load database table
dbtib<-read_tsv("BACT_DB_8.tsv", col_types = "ccdcdcdcddddddddffffffffffffffffffcddffffffdddddddddddddddddddd")
dbtib$sign.most<-factor(dbtib$sign.most, levels = c("pos","ns","neg"))
dbtib$len.cat<-factor(dbtib$len.cat, levels = c("[1-9]","[10-17]","[18-29]","[30-47]", "[48+]"))
dbtib

group.names<-c("neg"="NEG", "pos"="POS", "ns"="NS")

## Read gc content per position table
nuc.prob <- read_tsv(file = "NucleotideProbabilitiesRanseq-long.tsv", col_types = "cfd")
nuc.prob

## Load biased database simulation (same GC content/position as real DB)
largesim<- read_tsv("BACT_SimL_8.tsv")
largesim$len.cat<-factor(largesim$len.cat, levels = c("[1-9]","[10-17]","[18-29]","[30-47]", "[48+]"))
largesim

## Joint tables
jointib<- read_tsv("BAC_jointDB_8.tsv", col_types = "fccdcdcdcddddddddffffffffffffffffffcddffffffdddddddddddddddddddd")
jointib$len.cat<-factor(jointib$len.cat, levels = c("[1-9]","[10-17]","[18-29]","[30-47]", "[48+]"))
jointib$DB<-factor(jointib$DB, levels =c("Simulation", "Database"))
jointib

## Obtain peptide sets
peps.sim<-largesim[["pep.seq"]]
pepset.sim<-AAStringSet(x=peps.sim)
pepset.sim

peps<-dbtib[["pep.seq"]]
pepset<-AAStringSet(x=peps)
pepset

##Set figure output directory to no-mismatch directory
setwd("Figures-BL")

geofun<-function(x){
    nrow(dbtib)*(1-3/64)^(x-1)*3/64
}
##----------Fig2-LengthDistribution----------
pdf("Fig2-LengthDistribution.pdf", width = 8, height = 5)
ggplot(data = dbtib)+
    geom_histogram(mapping = aes(x=pep.len, fill="Database"), binwidth = 1, alpha=0.8)+
    geom_function(aes(colour="Black"), fun=geofun, linetype="solid", size=0.8)+
    coord_cartesian(xlim=c(0,65), ylim=c(0,800))+
    labs(x="Length (AA)", y="Number of sequences",caption=paste("n=", nrow(dbtib), "sequences"))+
    scale_fill_manual(values = colscale, labels=NULL, guide=NULL)+
    scale_colour_manual(values = colscale, labels="Geometric probability function", name=NULL)+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = c(0.1,-0.15), legend.background = element_blank())
dev.off()

##----------SupFig1-ClusterSizeDistribution----------
pdf("SupFig1-ClusterSizeDistribution.pdf", width = 8, height = 5)
ggplot(data = dbtib)+
    geom_histogram(mapping = aes(x=log10(cluster.size)), binwidth =0.01, size=0.5)+
    labs(x="Log10(Cluster Size)", y="Count",caption=paste("n=", nrow(dbtib), "sequences"))+
    theme_minimal(base_family = "Helvetica")
dev.off()

##----------SupFig2-GCcontent----------
## A. GC content in full seqs, ORF seqs and Random part of the seq (trimmed, full seqs)
gc<-ggplot(data = dbtib)+
    geom_density(aes(x=per.gc.rand, fill="RandomSeq"), linetype="blank", alpha=0.4)+
    geom_density(aes(x=per.gc.full, fill="Database"), linetype="blank", alpha=0.4)+
    geom_density(aes(x=per.gc.orf, fill="ORFs"), linetype="blank", alpha=0.3)+
    geom_vline(aes(xintercept=mean(per.gc.rand), colour="RandomSeq"),linetype="dashed", size=0.5)+
    geom_vline(aes(xintercept=mean(per.gc.full), colour="Database"),linetype="dashed", size=0.5)+
    geom_vline(aes(xintercept=mean(per.gc.orf), colour="ORFs"),linetype="dashed", size=0.5)+
    labs(x="GC%", y="Density", fill="",caption=paste("n=", nrow(dbtib)))+
    scale_colour_manual(values = colscale, name="Mean", labels=c("Full read", "ORF", "Random sequence"))+
    scale_fill_manual(values = colscale, name="Mean", labels=c("Full read", "ORF", "Random sequence"))+
    coord_cartesian(xlim=c(20,80))+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = c(0.3,0.8))
gc
# mean(dbtib$per.gc.full)
## 55.16352
# mean(dbtib$per.gc.orf)
## 53.04086
# mean(dbtib$per.gc.rand)
## 57.51811

## B. GC content per position
gc.pp<-ggplot(data=nuc.prob,mapping = aes(x=Position, y=Frequency, color=Nucleotide, group=Nucleotide))+
    geom_point(size=0.5, alpha=0.5)+
    geom_line()+
    geom_vline(aes(xintercept = 36), colour="darkgrey", size=0.5)+
    scale_color_manual(values=c("#89112d", "#000000", "#f7a707" ,"#1b69a4"), name=NULL)+
    labs(x="Position", y="Frequency",caption = paste0("n= ",sum(dbtib$rand.seq.len==150)))+
    scale_x_discrete(labels=c(rep("",35),"36",rep("",114)))+
    theme_minimal(base_family = "Helvetica")+
    theme(panel.grid.major.x = element_blank(), legend.position = "right", legend.text.align = 0, axis.ticks.length.x = unit(0,"mm"))
gc.pp

pdf("SupFig2-GCcontent.pdf", width = 8, height = 5)
ggarrange(gc, gc.pp, labels  =c("A", "B"), ncol = 2)
dev.off()

##----------SupFig-AAfrequencyDistributions----------
## Get amino acid frequencies from the joint database/simulation table
aa.freqs.join<-jointib%>%
    pivot_longer(W:P,names_to = "AA", values_to = "Frequency", values_drop_na = TRUE)%>%
    group_by(DB,AA)%>%
    mutate(Median=median(Frequency, na.rm = TRUE))%>%
    group_by(sign.most,AA)%>%
    mutate(MedianGroup=median(Frequency, na.rm = TRUE))%>%
    ungroup()
aa.freqs.join$AA<-factor(aa.freqs.join$AA, levels = myorder)
aa.freqs.join

## Get amino acid frequencies from the database table
aa.freqs.db<-dbtib%>%
    pivot_longer(W:P, names_to = "AA", values_to = "Frequency", values_drop_na = TRUE)%>%
    group_by(sign.most,AA)%>%
    mutate(Median=median(Frequency, na.rm = TRUE))
aa.freqs.db$AA<-factor(aa.freqs.db$AA, levels = myorder)
aa.freqs.db

## Get amino acid frequencies from the simulation table
aa.freqs.sim<-largesim%>%
    pivot_longer(W:P, names_to = "AA", values_to = "Frequency", values_drop_na = TRUE)%>%
    mutate(Median=median(Frequency, na.rm = TRUE))
aa.freqs.sim$AA<-factor(aa.freqs.sim$AA, levels = myorder)
aa.freqs.sim

## Compare the frequencies of the simulation and database for each amino acid
pdf(file = "SupFig-DensityPlotsaaFrequency.pdf")
ggplot(data=aa.freqs.join)+
    geom_density(aes(x=Frequency, colour=DB), size=0.8)+
    geom_vline(aes(xintercept=Median, colour=DB, linetype=DB), size=0.5)+
    labs(title = "Density")+
    theme_minimal()+
    scale_colour_manual(values = colscale, name=NULL)+
    scale_linetype_manual(values=c("43","73"), guide=NULL)+
    facet_wrap(facets = vars(AA), nrow = 5)

ggplot(data=aa.freqs.join,mapping = aes(x=Frequency, colour=DB))+
    stat_ecdf(geom = "step", pad=FALSE, size=1)+
    geom_vline(aes(xintercept=Median, colour=DB, linetype=DB), size=0.5)+
    labs(title = "Empirical cumulative distribution")+
    theme_minimal()+
    scale_colour_manual(values = colscale, name=NULL)+
    scale_linetype_manual(values=c("43","73"), guide=NULL)+
    facet_wrap(facets = vars(AA), nrow = 5)
dev.off()

##KS test to determine whether the Simulated distribution and the Empirical distribution of frequences are the same for each amino acid
## Calculate medians
medians<-NULL
## Calculate interquantile ranges
iqrs<-NULL
## Prepare result table
res<-NULL
for(i in myorder){
    print(i)
    x<-tapply(jointib[[i]], jointib$DB, median,na.rm=TRUE)
    medians<-bind_rows(medians, x)
    y<-tapply(jointib[[i]], jointib$DB, IQR, na.rm=TRUE)
    iqrs<-bind_rows(iqrs,y)
    db.aa<-dbtib%>%
        drop_na(i)
    sim.aa<-largesim%>%
        drop_na(i)
    a<-ks.test(db.aa[[i]], sim.aa[[i]])
    b<-tibble(Comparison="DB_vs_Simulation",AA=i,D=a$statistic, pvalue=a$p.value)
    res<-bind_rows(res,b)
}
res

## Correct for multiple testing
res<-bind_cols(res, p.adj=p.adjust(res$pvalue, method = "bonferroni"))
res
res<-res%>%
    mutate(label=case_when(p.adj<0.01 ~ "**", p.adj<0.05 ~ "*"))
res

## Check medians and interquartile ranges
medians<-bind_cols(AA=myorder, medians)
medians<-medians%>%
    mutate(.,diff=Simulation-Database)
medians
iqrs<-bind_cols(AA=myorder, iqrs)
iqrs<-iqrs%>%
    mutate(.,diff=Simulation-Database)
iqrs

##----------Fig3-AAfrequencies Plot AA distributions and expected values----------
## Get median values for the simulation table
aa_freqs.sim.med<-as_tibble(letterFrequency(pepset.sim,letters = myorder,as.prob = TRUE))
aa_freqs.sim.med<-aa_freqs.sim.med%>%
    na_if(0)
aa_freqs.sim.med<-apply(aa_freqs.sim.med, 2, median, na.rm=TRUE)
aa_freqs.sim.med<-tibble(AA=myorder, Frequency=aa_freqs.sim.med, label=res$label)
aa_freqs.sim.med

pdf("Fig3-AAfrequencies.pdf", width = 8, height = 5)
ggplot(aa.freqs.db, aes(x=AA, y=Frequency))+
    geom_violin(aes(fill="Low"), colour="lightgrey")+
    geom_boxplot(width=0.15, notch = TRUE, outlier.fill = "white", outlier.shape = 1, outlier.size = 0.7)+
    geom_errorbar(data=aa_freqs.sim.med, aes(ymax=Frequency, ymin=Frequency, colour="Black"), linetype="dashed")+
    geom_label(data = aa_freqs.sim.med, aes(label=label, y=0.55), size=6,label.size = 0, label.padding = unit(0, "lines"), hjust=.5, vjust=.5)+
    labs(x="Amino Acid", y="Frequency",caption = paste0("Total number of sequences=",nrow(dbtib),"\nBonferroni-Adjusted *P<0.05, **P<0.01"))+
    scale_colour_manual(values=colscale, name=NULL, label=c("Simulation median"))+
    scale_fill_manual(values=colscale, guide=NULL)+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position=c(0.1,-0.13), legend.spacing = unit(0, "cm"), plot.caption.position = "plot")


aa.freqs.join%>%
    select(DB,len.cat,AA,Frequency)%>%
    group_by(DB,AA,len.cat)%>%
    summarise(freq=mean(Frequency))%>%
    ggplot(aes(x=AA))+
    geom_col(aes(y=freq, fill=len.cat))+
    scale_fill_manual(values=colscalelen)+
    facet_wrap(facets = vars(DB))
    
aa.freqs.db%>%
    select(len.cat,sign.most,AA,Frequency)%>%
    group_by(AA,len.cat, sign.most)%>%
    summarise(freq=mean(Frequency))%>%
    ggplot(aes(x=AA))+
    geom_col(aes(y=freq, fill=len.cat), position="fill")+
    scale_fill_manual(values=colscalelen)+
    facet_wrap(facets = vars(sign.most))
dev.off()

##----------Intrinsic disorder (average)----------
## Get distributions for database and simulation
pdf("SupFig-IDSLen_DB_vs_Simulation.pdf", width = 8, height = 5)
ggplot(data = jointib)+
    geom_histogram(mapping = aes(x=iupred.short.avg, fill=len.cat),binwidth=0.01, position = "stack")+
    labs(x="Average IUPred score - short", y="Count" ,caption = paste0("n=",nrow(dbtib)))+
    scale_fill_manual(values = colscalelen, name="Peptide Length")+
    facet_grid(rows = vars(DB), scales = "free_y")+
    theme_minimal(base_family = "Helvetica")
dev.off()

##----------Plot IDS divided by length----------
##Other supplementary IDS data
pdf(file = "OtherSupplementaryDataIDS.pdf")
ggplot(data = dbtib)+
    geom_histogram(mapping = aes(x=iupred.short.avg),binwidth=0.01)+
    labs(x="Average IUPred score - short", y="Count" ,caption = paste0("Number of sequences=",nrow(dbtib)))+
    theme_minimal(base_family = "Helvetica")+
    facet_wrap(facets = vars(len.cat))

dbtib%>%
    filter(len.cat!="[1-9]")%>%
    ggplot()+
    geom_histogram(mapping = aes(x=iupred.short.avg, fill=len.cat),binwidth=0.01)+
    labs(x="Average IUPred score - short", y="Count" ,caption = paste0("Number of sequences=",nrow(dbtib)))+
    scale_fill_manual(values = colscalelen[-1], name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")

mean(dbtib[which(dbtib$len.cat!="[1-9]"),]$iupred.short.avg)
# 0.4494508

dbtib%>%
    group_by(len.cat)%>%
    summarise(mean=mean(iupred.short.avg), median=median(iupred.short.avg))

# # A tibble: 5 x 3
# len.cat  mean median
# * <fct>   <dbl>  <dbl>
# 1 [1-9]     0.947  0.988
# 2 10-17   0.677  0.707
# 3 18-29   0.479  0.472
# 4 30-47   0.362  0.347
# 5 [48+]     0.281  0.269

dbtib%>%
   filter(len.cat!="[1-9]")%>%
    group_by(sign.most)%>%
    summarise(mean=mean(iupred.short.avg), median=median(iupred.short.avg))

# A tibble: 4 x 3
# sign.most  mean median
# * <fct>     <dbl>  <dbl>
#     1 pos       0.502  0.493
# 2 ns        0.511  0.490
# 3 neg       0.343  0.305
# 4 NA        0.485  0.450

dbtib%>%
    group_by(sign.most)%>%
    summarise(mean=mean(iupred.short.avg), median=median(iupred.short.avg))
# # A tibble: 4 x 3
# sign.most  mean median
# * <fct>     <dbl>  <dbl>
# 1 pos       0.561  0.548
# 2 ns        0.611  0.603
# 3 neg       0.470  0.366
# 4 NA        0.578  0.551

##----------Plot IDS (fraction of disordered residues)----------
ggplot(data = dbtib)+
    geom_histogram(mapping = aes(x=iupred.short.frac, fill=len.cat),binwidth=0.02)+
    labs(x="Fraction of disordered resues IUPred - short", y="Count" ,caption = paste0("Number of sequences=",nrow(dbtib)))+
    scale_fill_manual(values = colscalelen, name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")

##----------Plot IDS (Average IDS IUpred long)----------
ggplot(data = dbtib)+
    geom_histogram(mapping = aes(x=iupred.long.avg, fill=len.cat),binwidth=0.01)+
    labs(title="Intrinsic disorder", x="Average IUPred score", y="Count" ,caption = paste0("n=",nrow(dbtib)))+
    scale_fill_manual(values = colscalelen, name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")

##----------Plot IDS (Average IDS IUpred long, divided by length category)----------
ggplot(data = dbtib)+
    geom_histogram(mapping = aes(x=iupred.long.avg),binwidth=0.01)+
    labs(x="Average IUPred score", y="Count" ,caption = paste0("Number of sequences=",nrow(dbtib)))+
    theme_minimal(base_family = "Helvetica")+
    facet_wrap(facets = vars(len.cat), scales = "free_y")

##----------Plot IDS (IUpred long, fraction of disordered residues)----------
ggplot(data = dbtib)+
    geom_histogram(mapping = aes(x=iupred.long.frac, fill=len.cat),binwidth=0.02)+
    labs(title="Intrinsic disorder", x="Fraction of disordered resues IUPred", y="Count" ,caption = paste0("Number of sequences=",nrow(dbtib)))+
    scale_fill_manual(values = colscalelen, name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")

dbtib%>%
    group_by(len.cat)%>%
    summarise(mean=mean(iupred.short.frac), median=median(iupred.short.frac))

# # A tibble: 5 x 3
# len.cat  mean median
# * <fct>   <dbl>  <dbl>
# 1 [1-9]     0.994  1    
# 2 10-17   0.795  1    
# 3 18-29   0.507  0.464
# 4 30-47   0.320  0.277
# 5 [48+]     0.207  0.185

##----------Plot IDS (IUpred long, fraction of disordered residues)----------
ggplot(data = dbtib)+
    geom_histogram(mapping = aes(x=iupred.long.frac),binwidth=0.01)+
    labs(x="Average IUPred score", y="Count" ,caption = paste0("Number of sequences=",nrow(dbtib)))+
    theme_minimal(base_family = "Helvetica")+
    facet_wrap(facets = vars(len.cat), scales = "free_y")
dev.off()

##----------Fig4-IDSvsPeptideLength----------
##Plot IDS distribution
ids<-ggplot(data = dbtib)+
    geom_histogram(mapping = aes(x=iupred.short.avg, fill=len.cat),binwidth=0.01, position = "stack")+
    labs(x="Average IUPred score", y="Count")+
    scale_fill_manual(values = colscalelen, name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = c(0.2,0.8))
ids

##Plot IDS vs length category
ids.pepcat<-ggplot(data = dbtib, mapping = aes(x=len.cat, y=iupred.short.avg))+
    geom_boxplot(aes(colour=len.cat, fill=len.cat), width=0.4)+
    labs(x="Peptide length category", y="Average IUPred score" )+
    scale_colour_manual(values = colscalelen, name="Peptide Length")+
    scale_fill_manual(values = colscalelen, name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position="none")+
    coord_flip()
ids.pepcat

##Get median coordinates from the plot build
dat<-ggplot_build(ids.pepcat)$data[[1]]
##Add median line as a white line in the middle
ids.pepcat<-ids.pepcat+ geom_segment(data=dat, aes(x=xmin, xend=xmax,y= middle, yend=middle), colour="white", size=1, alpha=0.6)
ids.pepcat

##Plot IDS vs peptide length
ids.peplen<-ggscatter(data = dbtib, x= "iupred.short.avg", y="pep.len", 
          add = "reg.line",
          conf.int = TRUE,
          color="len.cat",
          palette=colscalelen,
          size = 1,
          alpha=0.5)+
    stat_cor(aes(color = len.cat), label.x =0.70, label.y.npc = "top", method = "spearman")+
    labs(x="Average IUPred score", y="Predicted peptide length")+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = "none")
ids.peplen

##Plot IDS vs GC content
gc.ids<-ggplot(data = dbtib, aes(x= iupred.short.avg, y=per.gc.orf))+
    geom_point(aes(colour=len.cat))+
    geom_smooth(method = "lm", se = TRUE, aes(colour=len.cat))+
    stat_cor(aes(colour=len.cat), label.y.npc = "top", method = "spearman")+
    labs(x="Average IUPred score", y="GC content (ORF)" ,caption = paste0("n=",nrow(dbtib)))+
    scale_color_manual(values = colscalelen, guide=NULL)+
    theme_minimal(base_family = "Helvetica")
gc.ids

pdf("Fig4-IDS.pdf", width = 12, height = 10)
ggarrange(ids, ids.pepcat, ids.peplen, gc.ids, labels  =c("A", "B", "C", "D"), ncol = 2, nrow=2, align = "hv")
ggarrange(ids, ggarrange(ids.peplen, gc.ids, labels  =c("B", "C"), ncol = 2, align = "h"), nrow=2, labels = "A")
dev.off()

##----------Fig4-IDSvsPeptideLength (Alternative 2: Using fraction of disordered residues instead of average)----------
##4A.Plot IDS distribution
ids<-ggplot(data = dbtib)+
    geom_histogram(mapping = aes(x=iupred.short.frac, fill=len.cat),bins=40, position = "stack")+
    labs(x="Fraction of disordered residues", y="Count")+
    scale_fill_manual(values = colscalelen, name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = c(0.15,0.8))
ids

##Plot IDS vs length category
ids.pepcat<-ggplot(data = dbtib, mapping = aes(x=len.cat, y=iupred.short.frac))+
    geom_boxplot(aes(colour=len.cat, fill=len.cat), width=0.4)+
    labs(x="Peptide length category", y="Fraction of disordered residues" )+
    scale_colour_manual(values = colscalelen, name="Peptide Length")+
    scale_fill_manual(values = colscalelen, name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position="none")+
    coord_flip()
ids.pepcat
##Get median coordinates
dat<-ggplot_build(ids.pepcat)$data[[1]]
##Add median line as a white line in the middle
ids.pepcat<-ids.pepcat+ geom_segment(data=dat, aes(x=xmin, xend=xmax,y= middle, yend=middle), colour="white", size=1, alpha=0.6)
ids.pepcat

##Plot IDS vs peptide length
ids.peplen<-ggscatter(data = dbtib, x= "iupred.short.frac", y="pep.len", 
                      add = "reg.line",
                      conf.int = TRUE,
                      color="len.cat",
                      palette=colscalelen,
                      size = 1,
                      alpha=0.8,
                      position="jitter")+
    stat_cor(aes(color = len.cat), label.x =0.65, label.y.npc = 0.85, method = "spearman")+
    labs(x="Fraction of disordered residues", y="Predicted peptide length")+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = "none")
ids.peplen

##Plot IDS vs GC content
gc.ids<-ggplot(data = dbtib, aes(x= iupred.short.frac, y=per.gc.orf))+
    geom_point(aes(colour=len.cat))+
    geom_smooth(method = "lm", se = TRUE, aes(colour=len.cat))+
    stat_cor(aes(colour=len.cat),label.x.npc = 0.65, label.y.npc = 0.2, method = "spearman")+
    labs(x="Fraction of disordered residues", y="GC content (ORF)" ,caption = paste0("n=",nrow(dbtib)))+
    scale_color_manual(values = colscalelen, guide=NULL)+
    theme_minimal(base_family = "Helvetica")
gc.ids

pdf("Fig4-IDS(V2).pdf", width = 12, height = 10)
ggarrange(ids, ids.pepcat, ids.peplen, gc.ids, labels  =c("A", "B", "C", "D"), ncol = 2, nrow=2, align="hv")
ggarrange(ids, ggarrange(ids.peplen, gc.ids, labels  =c("B", "C"), ncol = 2, align = "h"), nrow=2, labels = "A")
dev.off()

##----------SupFig-NumberOfSequencesPerGroup----------
pdf("SupFig-NumberOfSequencesPerGroup.pdf", width = 8, height = 5)
ggplot(data = dbtib, aes(x=sign.most))+
    geom_bar(aes(fill=sign.most), width = 0.5)+
    scale_fill_manual(values = colscalegroup, na.value="slategray", guide=NULL)+
    geom_text(stat="count", aes(label=..count..), vjust=1.5, colour="white")+
    labs(x="General trend (>4 experiments)", y="Count")+
    scale_x_discrete(labels=group.names)+
    theme_minimal(base_family = "Helvetica")
dev.off()

##----------Fig5-LengthDistributionGroups----------
pdf("Fig5-LengthDistributionGroups.pdf", width = 9, height = 3)

annot<-tibble(text=c(paste("n=",sum(dbtib$sign.most=="pos",na.rm = TRUE)),paste("n=",sum(dbtib$sign.most=="ns",na.rm = TRUE)), paste("n=",sum(dbtib$sign.most=="neg",na.rm = TRUE))), x=c(55,55,55), y=c(650,650,650), sign.most=c("pos", "ns", "neg"))

dbtib.groups<-dbtib%>%
    filter(!is.na(sign.most))

dbtib.groups%>%
    ggplot(aes(x=pep.len))+
    geom_histogram( data=dbtib.groups %>% select(-sign.most), fill="lightgrey", binwidth = 1, na.rm = TRUE, alpha=0.4) +
    geom_histogram(mapping = aes(y=..count.., fill=sign.most),  binwidth = 1, na.rm = TRUE)+
    labs( x="Peptide Length (AA)", y="Number of sequences")+
   geom_density(aes(y=..density..*nrow(dbtib), colour=sign.most), linetype="dotted")+
    scale_x_continuous(breaks = c(0,20,40,60))+
    scale_fill_manual(values = colscalegroup, name=NULL, labels=group.names)+
    scale_color_manual(values = colscalegroup, guide=NULL)+
    geom_text(data = annot, aes(x=x, y=y,label=text), size=3)+
    coord_cartesian(xlim = c(4,65))+
    theme_minimal(base_family = "Helvetica")+
    facet_wrap(facets =vars(sign.most), labeller = as_labeller(group.names))+
    theme(legend.position = "none", panel.grid.minor = element_blank(), strip.background = element_blank(), panel.spacing.x = unit(5, "mm"), strip.text= element_text(face = "bold"))

dbtib.groups%>%
    ggplot(aes(x=pep.len))+
    geom_histogram( data=dbtib.groups %>% select(-sign.most), fill="lightgrey", binwidth = 1, na.rm = TRUE, alpha=0.4) +
    geom_histogram(mapping = aes(y=..count.., fill=sign.most),  binwidth = 1, na.rm = TRUE)+
    labs( x="Peptide Length (AA)", y="Number of sequences")+
    scale_x_continuous(breaks = c(0,20,40,60))+
    scale_fill_manual(values = colscalegroup, name=NULL, labels=group.names)+
    geom_text(data = annot, aes(x=x, y=y,label=text), size=3)+
    coord_cartesian(xlim = c(4,65))+
    theme_minimal(base_family = "Helvetica")+
    facet_wrap(facets =vars(sign.most), labeller = as_labeller(group.names))+
    theme(legend.position = "none", panel.grid.minor = element_blank(), strip.background = element_blank(), panel.spacing.x = unit(5, "mm"), strip.text= element_text(face = "bold"))

dbtib.groups%>%
    ggplot(aes(y=pep.len, x=sign.most))+
    geom_violin(mapping = aes(fill=sign.most, colour=sign.most), na.rm = TRUE, alpha=0.5)+
    geom_point(mapping = aes(fill=sign.most, colour=sign.most), na.rm = TRUE, alpha=0.5, position = position_jitterdodge())+
    labs(y="Peptide Length (AA)", x="Group")+
    scale_fill_manual(values = colscalegroup, name=NULL, labels=group.names)+
    scale_colour_manual(values = colscalegroup, name=NULL, labels=group.names)+
    theme_minimal(base_family = "Helvetica")

dev.off()

##----------Fig6-GCgroups----------
pdf("Fig6-GCgroups.pdf", width = 8, height = 5)
ggplot(data=dbtib)+
    geom_density(mapping = aes(x=per.gc.orf, fill=sign.most, colour=sign.most), alpha=0.2, size=1)+
    geom_vline(aes(xintercept=mean(per.gc.orf)), colour="black", linetype="dashed", size=0.5)+
    labs(x="GC% ORFs", y="Density",caption=paste("n=", nrow(dbtib)))+
    scale_fill_manual(values = colscalegroup, name=NULL, na.value=NA, na.translate=FALSE, labels=group.names)+
    scale_colour_manual(values = colscalegroup, name=NULL, na.value=NA, na.translate=FALSE, labels=group.names)+
    theme_minimal(base_family = "Helvetica")
dev.off()

##----------Fig7-OtherFeaturesGroups----------
##7A.Plot IDS distribution
ids.group<-dbtib%>%
    drop_na(sign.most)%>%
    ggplot()+
    geom_histogram(mapping = aes(x=iupred.short.avg, fill=len.cat),bins=80, position = "stack")+
    labs(x="Average IUPred2 score", y="Count")+
    scale_fill_manual(values = colscalelen, name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = "right")+
    facet_grid(rows = vars(sign.most), scales = "free_y", labeller = as_labeller(group.names))+
    theme(strip.background = element_rect(colour = "lightgray"), strip.text = element_text(face = "bold", size = 10))
ids.group

ids.group.ecdf<-dbtib%>%
    drop_na(sign.most)%>%
    ggplot(mapping = aes(x=iupred.short.avg, colour=len.cat))+
    stat_ecdf(geom = "step", pad=FALSE, size=0.5)+
    labs(title = "Empirical cumulative distribution", x="Average IUPred2 score", y="Count")+
    scale_colour_manual(values = colscalelen, name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")+
    theme(strip.background = element_rect(colour = "lightgray"), strip.text = element_text(face = "bold", size = 10), legend.position = "right", panel.grid.minor = element_blank())+
    facet_wrap(facets = vars(sign.most), nrow = 3, labeller = as_labeller(group.names), strip.position = "right")
ids.group.ecdf

ids.group2<-dbtib%>%
    drop_na(sign.most)%>%
    ggplot()+
    geom_histogram(mapping = aes(x=iupred.short.avg, fill=len.cat),bins=80, position = "stack")+
    labs(x="Average IUPred2 score", y="Count")+
    scale_fill_manual(values = colscalelen, name="Peptide Length")+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = "right")+
    facet_grid(cols = vars(sign.most), scales = "free_y", labeller = as_labeller(group.names))+
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10))
ids.group2

ids.group3<-ggplot(data = dbtib, mapping = aes(x=len.cat, y=iupred.short.avg))+
    geom_boxplot(aes(colour=sign.most, fill=sign.most), width=0.4, position = position_dodge(width = 0.8))+
    labs(x="Peptide length category", y="Average IUPred2 score" )+
    scale_colour_manual(values = colscalegroup, name="", na.translate=FALSE)+
    scale_fill_manual(values = colscalegroup, name="", na.translate=FALSE)+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position="right")+
    coord_flip()
ids.group3

##Get median coordinates from the plot build
dat<-ggplot_build(ids.group3)$data[[1]]
##Add median line as a white line in the middle
ids.group3<-ids.group3+ geom_segment(data=dat, aes(x=xmin, xend=xmax,y= middle, yend=middle), colour="white", size=1, alpha=0.6)
ids.group3

dbtib%>%
    drop_na(sign.most)%>%
    ggplot()+
    geom_density(aes(x=iupred.short.avg, y=..scaled.., color=sign.most), size=1)+
    labs( x="Average IUPred2 score", y="Density", caption=paste("n=", nrow(dbtib)))+
    scale_color_manual(values = colscalegroup, labels=group.names, name=NULL)+
    theme_minimal(base_family = "Helvetica")

## Table of p values
## KS test for comparison between group distributions
res2<-NULL

for(i in myorder){
    print(i)
    po<-dbtib[which(dbtib$sign.most=="pos"),]%>%
        drop_na(i)
    ne<-dbtib[which(dbtib$sign.most=="neg"),]%>%
        drop_na(i)
    ns<-dbtib[which(dbtib$sign.most=="ns"),]%>%
        drop_na(i)
    pne<-ks.test(po[[i]], ne[[i]])
    pns<-ks.test(po[[i]], ns[[i]])
    nne<-ks.test(ne[[i]], ns[[i]])
    b<-tibble(Comparison=c("Up_vs_Down","Up_vs_NS","Down_vs_NS"), AA=i,D=c(pne$statistic, pns$statistic, nne$statistic), pvalue=c(pne$p.value, pns$p.value,nne$p.value))
    res2<-bind_rows(res2,b)
}
res2
res2<-bind_cols(res2, p.adj=p.adjust(res2$pvalue, method = "bonferroni"))
res2<-res2%>%
    mutate(label=case_when(p.adj<0.01 ~ "**", p.adj<0.05 ~ "*"))
res2
res2$label
res3<-bind_rows(res,res2)%>%
    pivot_wider(id_cols = AA, names_from = Comparison, values_from = D:label)%>%
    select(AA, ends_with("DB_vs_Simulation"), ends_with("Up_vs_Down"), ends_with("Down_vs_NS"), ends_with("Up_vs_NS"))
res3

write_tsv(res3, file="pvaluesAAfreq.tsv")

aa.group<-ggplot(aa.freqs.db[which(!is.na(aa.freqs.db$sign.most)),], aes(x=AA, y=Frequency))+
    geom_violin(aes(fill=sign.most), colour=NA, position = position_dodge(width=0.8), alpha=0.5)+
    geom_boxplot(aes(colour=sign.most), width=0.15,position = position_dodge(0.8), outlier.shape = 1, outlier.size = 0.7)+
    geom_errorbar(data=aa_freqs.sim.med, aes(ymax=Frequency, ymin=Frequency), colour="black", linetype="dashed")+
    geom_label(data = res2, aes(label=label, y=c(0.55)), size=8,label.size = 0, label.padding = unit(0, "lines"), hjust=.5, vjust=.5)+
    labs(x="Amino Acid", y="Frequency",caption = paste0("n=",nrow(dbtib[which(!is.na(dbtib$sign.most)),]),"\nBonferroni-Adjusted **P<0.01"))+
    scale_colour_manual(values=colscalegroup, name=NULL, label=group.names)+
    scale_fill_manual(values=colscalegroup, name=NULL, label=group.names)+
    coord_cartesian(ylim = c(0,0.55))+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = c(0.5,-0.1), legend.direction = "horizontal")
aa.group

sum.aa<-aa.freqs.db%>%
    group_by(sign.most, AA)%>%
    drop_na(sign.most)%>%
    summarise(mean=mean(Frequency, na.rm=TRUE), sd=sd(Frequency, na.rm=TRUE))
sum.aa
aa.group2<-ggplot(sum.aa, aes(x=AA, y=mean,fill=sign.most))+
    geom_col(colour=NA, position = position_dodge(width=0.8))+
    geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd), colour="black", position = position_dodge(width=0.8))+
    labs(x="Amino Acid", y="Frequency",caption = paste0("n=",nrow(dbtib[which(!is.na(dbtib$sign.most)),])))+
    scale_fill_manual(values=colscalegroup, name=NULL, label=group.names)+
    theme_minimal(base_family = "Helvetica")

aa.group2

## Compare the amino acid frequencies of each amino acid for the groups
aa.dens.group<-ggplot(data=aa.freqs.join[which(!is.na(aa.freqs.join$sign.most)),])+
    geom_line(aes(x=Frequency, colour=sign.most), stat="density", size=0.8, alpha=0.7)+
    geom_vline(aes(xintercept=MedianGroup, colour=sign.most, linetype=sign.most), size=0.5)+
    labs(title = "Density")+
    theme_minimal()+
    scale_colour_manual(values = colscalegroup, name=NULL, labels=group.names)+
    scale_linetype_manual(values=c("43","53","63"), guide=NULL)+
    facet_wrap(facets = vars(AA), nrow = 5)
aa.dens.group

aa.ecdf.group<-ggplot(data=aa.freqs.join[which(!is.na(aa.freqs.join$sign.most)),],mapping = aes(x=Frequency, colour=sign.most))+
    stat_ecdf(geom = "step", pad=FALSE, size=1, alpha=0.7)+
    geom_vline(aes(xintercept=MedianGroup, colour=sign.most, linetype=sign.most), size=0.5)+
    labs(title = "Empirical cumulative distribution")+
    theme_minimal()+
    scale_colour_manual(values = colscalegroup, name=NULL)+
    scale_linetype_manual(values=c("43","53","63"), guide=NULL)+
    facet_wrap(facets = vars(AA), nrow = 5)
aa.ecdf.group

pdf("Fig7-GroupFeatures.pdf", width = 10, height = 6)
ids.group
ids.group.ecdf
ids.group2
ids.group3
aa.group
aa.dens.group
aa.ecdf.group


ggarrange(ids.group,aa.group, nrow=2, heights = c(1,2), labels = c("A", "B"))
dev.off()

##--------Fig8.AggregationGroups--------
## 
aggr.length<-ggplot(data = dbtib)+
    geom_density(mapping = aes(x=pasta.best.e, colour=len.cat),na.rm = TRUE, size=1)+
    geom_vline(xintercept = -5, linetype="dashed")+
    labs(x="Best pairing aggregation energy (PEU)", y="Density",caption=paste("n=", nrow(dbtib)))+
    scale_colour_manual(values = colscalelen, name=NULL, labels=group.names)+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = c(0.2,0.7), plot.caption = element_text(hjust = 0))
aggr.length

aggr.length2<-ggplot(data = dbtib)+
    geom_density(mapping = aes(x=pasta.best.e, colour=sign.most),na.rm = TRUE, size=0.6)+
    geom_vline(xintercept = -5, linetype="dashed")+
    labs(x="Best pairing aggregation energy (PEU)", y="Density",caption=paste("n=", nrow(dbtib[which(!is.na(dbtib$sign.most)),])))+
    scale_colour_manual(values = colscalegroup, name=NULL, na.translate=FALSE, labels=group.names)+
    theme_minimal(base_family = "Helvetica")+
    facet_wrap(facets = vars(len.cat), nrow = 5)+
    theme(plot.caption = element_text(hjust = 0))
aggr.length2

aggr.length3<-dbtib%>%
    filter(!is.na(sign.most))%>%
    ggplot()+
    geom_density(mapping = aes(x=pasta.best.e, colour=sign.most), size=1)+
    geom_vline(xintercept = -5, linetype="dashed")+
    labs(x="Best pairing aggregation energy (PEU)", y="Density",caption=paste("n=", nrow(dbtib[which(!is.na(dbtib$sign.most)),])))+
    scale_colour_manual(values = colscalegroup, name=NULL, labels=group.names)+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = c(0.2,0.7), plot.caption = element_text(hjust = 0))
aggr.length3

aggr.length3.2<-dbtib%>%
    filter(!is.na(sign.most))%>%
    ggplot()+
    geom_histogram(mapping = aes(x=pasta.best.e, fill=sign.most), binwidth = 0.08, position = position_dodge())+
    geom_vline(xintercept = -5, linetype="dashed")+
    labs(x="Best pairing aggregation energy (PEU)", y="Density",caption=paste("n=", nrow(dbtib[which(!is.na(dbtib$sign.most)),])))+
    scale_fill_manual(values = colscalegroup, name=NULL, labels=group.names)+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = c(0.2,0.7), plot.caption = element_text(hjust = 0))+
    facet_grid(facets = vars(sign.most), scales = "free_y")
aggr.length3.2

aggr.length4<-dbtib%>%
    filter(!is.na(sign.most))%>%
    ggplot()+
    geom_density(mapping = aes(x=pasta.best.e, colour=len.cat),na.rm = TRUE, size=0.6)+
    geom_vline(xintercept = -5, linetype="dashed")+
    labs(x="Best pairing aggregation energy (PEU)", y="Density",caption=paste("n=", nrow(dbtib[which(!is.na(dbtib$sign.most)),])))+
    scale_colour_manual(values = colscalelen, name=NULL, na.translate=FALSE, labels=group.names)+
    theme_minimal(base_family = "Helvetica")+
    facet_wrap(facets = vars(sign.most), nrow = 3, labeller = as_labeller(group.names) )+
    theme(plot.caption = element_text(hjust = 0))
aggr.length4

agg.pepcat<-dbtib%>%
    filter(!is.na(sign.most))%>%
    ggplot(aes(x=len.cat, y=pasta.best.e))+
    geom_boxplot(aes(colour=sign.most, fill=NA), width=0.3, position = position_dodge(width = 0.5))+
    labs(x="Peptide length category", y="Best pairing aggregation energy (PEU)" )+
    geom_hline(yintercept = -5, linetype="dashed")+
    scale_colour_manual(values = colscalegroup, name="Peptide Length")+
    scale_fill_manual(values = colscalegroup, name="Group")+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position="none")+
    coord_flip()
agg.pepcat

dbtib%>%
    filter(!is.na(sign.most))%>%
    ggplot(aes(x=sign.most, y=pasta.best.e))+
    geom_boxplot(aes(colour=len.cat), width=0.3, position = position_dodge(width = 0.4))+
    labs(x="Peptide length category", y="Best pairing aggregation energy (PEU)" )+
    geom_hline(yintercept = -5, linetype="dashed")+
    scale_colour_manual(values = colscalelen )+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position="none")

pdf("Fig8-AggregationGroups.pdf", width = 10, height = 5)
ggarrange(aggr.length3,agg.pepcat, ncol=2, widths=c(1,1) ,heights = c(2,2), labels = c("A", "B"))
ggarrange(aggr.length3,aggr.length4, ncol=2, widths=c(2,1) ,heights = c(3,3), labels = c("A", "B"))
ggarrange(aggr.length,aggr.length2, ncol=2, widths=c(2,1) ,heights = c(3,3), labels = c("A", "B"))
aggr.length
aggr.length2
aggr.length3
aggr.length4
agg.pepcat
dev.off()

dbtib%>%
    group_by(len.cat)%>%
    summarise(mean=mean(pasta.best.e), median=median(pasta.best.e))
# # A tibble: 5 x 3
# len.cat   mean median
# * <fct>    <dbl>  <dbl>
# 1 [1-9]     -0.273  0.121
# 2 10-17   -2.07  -1.77 
# 3 18-29   -3.38  -3.12 
# 4 30-47   -4.41  -4.15 
# 5 [48+]     -5.63  -5.21 

dbtib%>%
    group_by(sign.most)%>%
    summarise(mean=mean(pasta.best.e), median=median(pasta.best.e))
# # A tibble: 4 x 3
# sign.most  mean median
# * <fct>     <dbl>  <dbl>
# 1 pos       -2.76  -2.73
# 2 ns        -2.58  -2.31
# 3 neg       -3.96  -4.07
# 4 NA        -2.96  -2.80

sum(dbtib$pep.len<6&dbtib$pasta.best.e<=-2)
# 0
table(dbtib[which(dbtib$pasta.best.e>=0),"len.cat"])
# [1-9] 10-17 18-29 30-47   [48+] 
# 635   135    25     0     0 

