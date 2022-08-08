## Created by: Johana Fajardo (09.2021)
## Last modified: (29.09.2021)
## Changes:
## Created
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
colscale0<-c("#f7a707","#f5d36c","#f5e8c1","#3ecad7","#2285d0","#2d52d8","#57359c","#0c324f","#000000","#bcbdc1")
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

## Load biased database simulation (same GC content/position as real DB)
largesim<- read_tsv("BACT_SimL_8.tsv")
largesim$len.cat<-factor(largesim$len.cat, levels = c("[1-9]","[10-17]","[18-29]","[30-47]", "[48+]"))
largesim

## Count how many sequences will be removed from the require by this operation
sum(dbtib$pep.len<=4)
## 200
sum(largesim$pep.len<=4)
## 36657

dbtib<-dbtib[which(dbtib$pep.len>4),]
dbtib

largesim<-largesim[which(largesim$pep.len>4),]
largesim

## Obtain peptide sets and aa frequencies for only the random part of the sequences
## For simulation
peps.sim<-largesim[["pep.seq"]]
mywidths1<-NULL
for(i in 1:length(peps.sim)){
    if(width(peps.sim[i])>54){
        mywidths1[i]<-50
    }else{
        mywidths1[i]<-width(peps.sim[i])-4
    }
}
mywidths1
pepset.sim<-AAStringSet(x=peps.sim,start = 5, width = mywidths1)
pepset.sim

j.peps.sim<-paste0(unlist(pepset.sim,recursive = TRUE), collapse = "")
j.pepset.sim<-AAStringSet(j.peps.sim)
sim.aa.freqs<-as_tibble(letterFrequency(j.pepset.sim, letters = myorder, as.prob = TRUE))%>%
    bind_cols(Group="Simulation",.)

freq.aa2<-as_tibble(letterFrequency(pepset.sim,letters = myorder,as.prob = TRUE))%>%
    na_if(0)
largesim<-bind_cols(largesim[1:17],freq.aa2)
largesim
hist(width(pepset.sim), breaks=50)

## For database
peps<-dbtib[["pep.seq"]]
mywidths<-NULL
for(i in 1:length(peps)){
    if(width(peps[i])>54){
        mywidths[i]<-50
    }else{
        mywidths[i]<-width(peps[i])-4
    }
}
mywidths
pepset<-AAStringSet(x=peps,start = 5, width = mywidths)
pepset
hist(width(pepset), breaks = 50)

j.peps<-paste0(unlist(pepset,recursive = TRUE), collapse = "")
j.pepset<-AAStringSet(j.peps)
letterFrequency(j.pepset, letters = myorder, as.prob = TRUE)
db.aa.freqs<-as_tibble(letterFrequency(j.pepset, letters = myorder, as.prob = TRUE))%>%
    bind_cols(Group="Database",.)

posind<-which(dbtib$sign.most=="pos")
negind<-which(dbtib$sign.most=="neg")
nsind<-which(dbtib$sign.most=="ns")

j.peps.pos<-paste0(unlist(pepset[posind],recursive = TRUE), collapse = "")
j.pepset.pos<-AAStringSet(j.peps.pos)
letterFrequency(j.pepset.pos, letters = myorder, as.prob = TRUE)
db.aa.freqs.pos<-as_tibble(letterFrequency(j.pepset.pos, letters = myorder, as.prob = TRUE))%>%
    bind_cols(Group="POS",.)

j.peps.neg<-paste0(unlist(pepset[negind],recursive = TRUE), collapse = "")
j.pepset.neg<-AAStringSet(j.peps.neg)
letterFrequency(j.pepset.neg, letters = myorder, as.prob = TRUE)
db.aa.freqs.neg<-as_tibble(letterFrequency(j.pepset.neg, letters = myorder, as.prob = TRUE))%>%
    bind_cols(Group="NEG",.)

j.peps.ns<-paste0(unlist(pepset[nsind],recursive = TRUE), collapse = "")
j.pepset.ns<-AAStringSet(j.peps.ns)
letterFrequency(j.pepset.ns, letters = myorder, as.prob = TRUE)
db.aa.freqs.ns<-as_tibble(letterFrequency(j.pepset.ns, letters = myorder, as.prob = TRUE))%>%
    bind_cols(Group="NS",.)

aa.freqs.all<-bind_rows(sim.aa.freqs, db.aa.freqs, db.aa.freqs.pos, db.aa.freqs.neg, db.aa.freqs.ns)%>%
    pivot_longer(names_to = "AA", cols = W:P)%>%
    pivot_wider(names_from = Group)
aa.freqs.all

aa.perc.all<-aa.freqs.all%>%
    mutate(across(Simulation:NS, ~ .x*100))

aa.perc.vsdb<-aa.perc.all%>%
    mutate(POS=POS-Database, NEG=NEG-Database, NS=NS-Database)%>%
    select(AA, POS, NEG, NS)%>%
    pivot_longer(cols = POS:NS, names_to = "Comparison", values_to = "Difference")%>%
    mutate(Sign=sign(Difference))
aa.perc.vsdb$Comparison<-factor(aa.perc.vsdb$Comparison, levels=c("POS", "NS", "NEG"))
aa.perc.vsdb$AA<-factor(aa.perc.vsdb$AA, levels=myorder)
aa.perc.vsdb$Sign<-factor(aa.perc.vsdb$Sign)
aa.perc.vsdb

aa.perc.vsgroups<-aa.perc.all%>%
    mutate(POS_minus_NEG=POS-NEG, NEG_minus_NS=NEG-NS, POS_minus_NS=POS-NS)%>%
    select(AA, contains("minus"))%>%
    pivot_longer(cols = contains("minus"), names_to = "Comparison", values_to = "Difference")%>%
    mutate(Sign=sign(Difference))
aa.perc.vsgroups$AA<-factor(aa.perc.vsdb$AA, levels=myorder)
aa.perc.vsgroups$Comparison<-factor(aa.perc.vsgroups$Comparison, levels=c("POS_minus_NEG", "POS_minus_NS", "NEG_minus_NS"))
aa.perc.vsgroups$Sign<-factor(aa.perc.vsgroups$Sign)
aa.perc.vsgroups

x<-aa.perc.vsdb%>%
    ggplot(aes(x=AA))+
    geom_col(aes(y=Difference, fill=Comparison), position = position_dodge())+
    geom_vline(xintercept=seq(0,20)+.5,color="lightgray", size=0.2)+
    labs(y="% Difference\n(Group-Database)")+
    coord_cartesian(ylim = c(-1.3, 1.3))+
    theme_minimal()+
    scale_fill_manual(values = colscalegroup, name=NULL)+
    theme(legend.key.size = unit(0.35, "cm"), legend.position = c(0,0), legend.justification = c(0,0), legend.background = element_rect(fill="white", colour = NA), panel.grid.minor = element_blank(), panel.grid.major=element_blank())
x

y<-aa.perc.vsgroups%>%
    ggplot(aes(x=AA))+
    geom_col(aes(y=Difference, fill=Comparison), position = position_dodge())+
    geom_vline(xintercept=seq(0,20)+.5,color="lightgray", size=0.2)+
    labs(y="% Difference")+
    coord_cartesian(ylim = c(-1.3, 1.3))+
    theme_minimal()+
    scale_fill_manual(values = c("#7BA807", "#F5DE98", "#389FA8"), name=NULL)+
    theme(legend.key.size = unit(0.35, "cm"), legend.position = c(0,0), legend.justification = c(0,0), legend.background = element_rect(fill="white", colour = NA), panel.grid.minor = element_blank(), panel.grid.major=element_blank())
y

v<-aa.perc.vsgroups%>%
    ggplot(aes(x=AA, y=abs(Difference)))+
    geom_point(aes(colour=Comparison, shape=Sign), size=4, stroke=2, alpha=0.85)+
    scale_color_manual(values =c("#7BA807", "#F5DE98", "#389FA8"), name=NULL)+
    guides(shape=guide_legend(reverse=TRUE,override.aes = list(size=c(1.5,1), stroke=1)), colour=guide_legend(override.aes = list(size=2, stroke=1)))+
    scale_shape_manual(values = c("triangle down filled", "triangle"), name=NULL, labels=c("Decreased", "Increased"))+
    labs(y="Absolute difference")+
    theme_minimal()+
    theme(panel.grid.minor = element_blank(), legend.key.size = unit(1, "pt"),legend.position = c(0,1), legend.justification = c(0,1), legend.spacing = unit(0,"cm"), legend.background = element_rect(fill="white", colour = NA))+
    coord_cartesian(ylim = c(0,1.3))
v

w<-aa.perc.vsdb%>%
    ggplot(aes(x=AA, y=abs(Difference)))+
    geom_point(aes(colour=Comparison, shape=Sign), size=4, stroke=2, alpha=0.85)+
    scale_shape_manual(values = c("triangle down filled", "triangle"), labels= c("Decreased", "Increased"), name=NULL)+
    geom_vline(xintercept=seq(0,20)+.5,color="lightgray", size=0.2)+
    scale_color_manual(values = colscalegroup, name=NULL)+
    labs(y="Absolute difference\n(Group-Database)")+
    guides(shape=guide_legend(order=2, reverse=TRUE, override.aes = list(size=c(1.5,1), stroke=1)), colour=guide_legend(order=1, override.aes = list(size=2, stroke=1)))+
    theme_minimal()+
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),legend.key.size = unit(1, "pt"), legend.position = c(0,1), legend.justification = c(0,1), legend.background = element_rect(fill="white", colour=NA), legend.spacing = unit(0, "cm"))+
    coord_cartesian(ylim = c(0,1.3))
w

pdf("AAcomparisons.pdf", width = 8, height = 4)
x
y
dev.off()

pdf("AAcomparisons_Combined.pdf", width = 8, height = 8)
ggarrange(x, y, nrow = 2, labels = "AUTO", align = "v")
dev.off()

pdf("AAcomparisons_Absolute_Combined.pdf", width = 8, height = 8)
ggarrange(w, v, nrow = 2, labels = "AUTO", align = "v")
dev.off()

pdf("AAcomparisons_Absolute.pdf", width = 8, height = 4)
w
v
dev.off()

write_tsv(aa.freqs.all, file="AA-Freqs-all.tsv")

freq.aa<-as_tibble(letterFrequency(pepset,letters = myorder,as.prob = TRUE))%>%
    na_if(0)
dbtib<-bind_cols(dbtib[1:43],freq.aa)
dbtib

## Join tables
dbtib<-bind_cols(DB="Database", dbtib)
dbtib

largesim<-bind_cols(DB="Simulation",largesim)
largesim

jointib<-bind_rows(dbtib,largesim)%>%
    mutate(.,across(.cols = c(DB, len.cat), as.factor))
jointib$DB<-factor(jointib$DB, levels = c("Simulation", "Database"))
jointib

##Set figure output directory to no-mismatch directory
setwd("Figures-BL")

geofun<-function(x){
    nrow(dbtib)*(1-3/64)^(x-1)*3/64
}


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
## Get median values for the simulaation table
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
dev.off()


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
    geom_boxplot(aes(colour=sign.most), width=0.15,position = position_dodge(0.8), outlier.shape = 1, outlier.size = 0.7, size=0.35)+
    geom_errorbar(data=aa_freqs.sim.med, aes(ymax=Frequency, ymin=Frequency), colour="black", linetype="dashed")+
 #   geom_label(data = res2, aes(label=label, y=c(0.55)), size=8,label.size = 0, label.padding = unit(0, "lines"), hjust=.5, vjust=.5)+
    labs(x="Amino Acid", y="Frequency",caption = paste0("n=",nrow(dbtib[which(!is.na(dbtib$sign.most)),])))+
    scale_colour_manual(values=colscalegroup, name=NULL, label=group.names)+
    scale_fill_manual(values=colscalegroup, name=NULL, label=group.names)+
    coord_cartesian(ylim = c(0,0.55))+
    theme_minimal(base_family = "Helvetica")+
    theme(legend.position = c(0.5,-0.2), legend.direction = "horizontal")
aa.group

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

pdf("Fig7-GroupFeatures.pdf", width = 10, height = 5)
aa.group
aa.dens.group
aa.ecdf.group
dev.off()

a<-dbtib%>%
    group_by(sign.most)%>%
    summarise(across(.cols = W:P, median, na.rm=TRUE))%>%
    rename(group=sign.most)

b<-jointib%>%
    group_by(DB)%>%
    summarise(across(.cols = W:P, median, na.rm=TRUE))%>%
    rename(group=DB)

c<-bind_rows(a,b)%>%
    pivot_longer(cols = W:P)%>%
    pivot_wider(names_from = group)


write_tsv(c, "AA-frequency-medians.tsv")
