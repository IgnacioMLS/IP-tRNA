
list.of.packages <- c("GenomicAlignments","GenomicFeatures", "Rsamtools","edgeR",
                      "magick", "Glimma","reshape2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  install.packages(new.packages, repos = "http://cran.us.r-project.org")
}


suppressPackageStartupMessages({
library("Rsubread")
library( "GenomicAlignments" )
library( "GenomicFeatures" )
library( "Rsamtools" )
library("edgeR")
library("vsn")
library("magick")
library("Glimma")
library("pheatmap")
library('tidyr')
library('stringr')
library(gtools)
library(stringr)
library(gdata)
library(dplyr)
library(plyr)
library(data.table)
library(heatmaply)
library(gridExtra)
library(ggplot2)
})


sample_data = read.table("../Fastq_downloaded/sample_data.txt", header=T,
                         stringsAsFactors = T)

groups = levels(sample_data$Condition)
dir = "../Results/Modification_ratio_plots"

aas = c("Ala", "Arg", "Asn", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "iMet", 
        "Leu", "Lys", "Met", "Phe", "Pro", "SeC", "Ser", "Thr", "Trp", "Tyr", "Val")
files <- list.files(path = dir, pattern = ".txt")
if (!file.exists("../Results/Modification_ratio_plots/Comparison")){
  dir.create("../Results/Modification_ratio_plots/Comparison", recursive=TRUE)
}

number_positions = 4
aa_info = data.frame(pos="", gene="",pvalue="")
for(aa in aas){
  files_aa = grep(pattern=aa, files, value=T)
  genes_aa = levels(as.factor(grep(pattern=groups[1], files_aa, value=T)))
  genes_aa2 = levels(as.factor(grep(pattern=groups[2], files_aa, value=T)))
  
  gene_levels = substr(genes_aa, start=nchar(groups[1])+2, nchar(genes_aa))
  gene_levels2 = substr(genes_aa2, start=nchar(groups[2])+2, nchar(genes_aa2))
  # Only use the genes that have data in both groups.
  gene_levels = gene_levels[gene_levels %in% gene_levels2]
  
  positions = c("9","27","32", "34", "37","47", "58")
  position_info = matrix(ncol=length(positions), nrow=1)
  colnames(position_info) = positions
  rownames(position_info) = "pvalue"
  real_genes= substring(gene_levels, 1, nchar(gene_levels)-4)
  na =rep(NA,length(gene_levels))
  #
  nas = na
  for(i in 1:(length(positions)-1)){
    nas = cbind(nas,na)
  }
  gene_info = t(nas)
  colnames(gene_info) = real_genes
  rownames(gene_info) = positions
  cont = 0
  
  for(gene in gene_levels){
    cont = cont +1
    data_group1 = read.table(file=paste0(dir,"/",groups[1],"_",gene), header=T)
    data_group2 = read.table(file=paste0(dir,"/",groups[2],"_",gene), header=T)
    colnames(data_group1)[2] = "Corrected_pos" 
    colnames(data_group2)[2] = "Corrected_pos" 
    
    data_group1 = data.frame(data_group1, Group = groups[1])
    data_group2 = data.frame(data_group2, Group =groups[2])
    data_tot = rbind(data_group1, data_group2)
    real_gene = substring(gene,1,nchar(gene)-4)
    data_tot$Modification_ratio
    p = ggplot(data = data_tot, aes(x=Corrected_pos, y=as.numeric(Modification_ratio), 
                                color=Group)) + geom_line()+
      labs(title='Modification ratio (Relative to gene coverage)', subtitle=real_gene,
           x="Position", y ="Modification ratio (Relative to gene coverage)")+
      scale_fill_manual(values= c("#0000FF", "#FF0000"),
                        breaks=c(groups[1], groups[2]))+
      scale_x_discrete(limits = data_group1$Corrected_pos)+
      theme(axis.text=element_text(size=16, color="black"),
            axis.text.x = element_text(size=12,angle=90),
            axis.title= element_text(size=20), legend.text = element_text(size = 16),
            legend.title = element_text(size=20), plot.subtitle = element_text(size=19),
            plot.title =element_text(size=20))
    
    file_name = paste0("Modification_ratio_by_pos_", real_gene, ".jpeg")
    ggsave(p, file=file_name, width = 35, 
           height = 25, units = "cm", path="../Results/Modification_ratio_plots/Comparison" )
  
    ###FISHER TEST ####
    for(i in 1:length(positions)){
      group1_pos<- as.data.frame(data_group1[data_group1$Corrected_pos == positions[i], ])
      group2_pos<- as.data.frame(data_group2[data_group2$Corrected_pos == positions[i], ])
      position_pos_test <- matrix(c(group1_pos$Base_coverage,group1_pos$Mod_bases,
                                   group2_pos$Base_coverage,group2_pos$Mod_bases),
                                 ncol=2,byrow=TRUE)
      
      rownames(position_pos_test)<-groups
      colnames(position_pos_test)<-c("Coverage","Modification")
      test_pos <- fisher.test(position_pos_test)
      p_pos <- test_pos$p.value
      position_info[i] = p_pos
    }
    gene_info[,cont] = position_info
  }
  indx = gene_info<0.05
  pos = rownames(gene_info)[row(gene_info)*indx]
  gene = colnames(gene_info)[col(gene_info)*indx]
  pvalue = gene_info[indx]
  gene_info_final = data.frame(pos, gene, pvalue)
  aa_info = rbind(aa_info, gene_info_final)
  }

aa_info = aa_info[aa_info$pos!="",]
write.table(aa_info, file = "../Results/R_files/Modification_test.txt", 
            quote=FALSE, row.names = F)

