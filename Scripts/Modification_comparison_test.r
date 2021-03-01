
list.of.packages <- c("GenomicAlignments","GenomicFeatures", "Rsamtools","edgeR",
                      "magick", "Glimma","reshape2", "RColorBrewer")
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
  library('RColorBrewer')
})


sample_data = read.table("../Fastq_downloaded/sample_data.txt", header=T,
                         stringsAsFactors = T)

groups = levels(sample_data$Condition)
dir = "../Results/Modification_ratio_plots"
dir_scripts = getwd()

aas = c("Ala", "Arg", "Asn", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "iMet", 
        "Leu", "Lys", "Met", "Phe", "Pro", "SeC", "Ser", "Thr", "Trp", "Tyr", "Val")
files <- list.files(path = dir, pattern = ".txt")
if (!file.exists("../Results/Modification_ratio_plots/Comparison")){
  dir.create("../Results/Modification_ratio_plots/Comparison", recursive=TRUE)
}

aa_info = data.frame(pos="", gene="",pvalue="")
number_genes = 1
for(aa in aas){
  files_aa = grep(pattern=aa, files, value=T)
  genes_aa = levels(as.factor(grep(pattern=groups[1], files_aa, value=T)))
  genes_aa2 = levels(as.factor(grep(pattern=groups[2], files_aa, value=T)))
  
  gene_levels = substr(genes_aa, start=nchar(groups[1])+2, nchar(genes_aa))
  gene_levels2 = substr(genes_aa2, start=nchar(groups[2])+2, nchar(genes_aa2))
  # Only use the genes that have data in both groups.
  gene_levels = gene_levels[gene_levels %in% gene_levels2]
  
  positions = c("9","26","32", "34", "37","47", "58")
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
  
  for(gene in real_genes){
    cont = cont +1
    data_group1 = read.table(file=paste0(dir,"/",groups[1],"_",gene, ".txt"), header=T)
    data_group2 = read.table(file=paste0(dir,"/",groups[2],"_",gene, ".txt"), header=T)
    colnames(data_group1)[2] = "Corrected_pos" 
    colnames(data_group2)[2] = "Corrected_pos" 
    
    data_group1 = data.frame(data_group1, Group = groups[1])
    data_group2 = data.frame(data_group2, Group =groups[2])
    data_tot = rbind(data_group1, data_group2)
    data_tot$Group = factor(data_tot$Group)
    p = ggplot(data = data_tot, aes(x=Corrected_pos, y=as.numeric(Modification_ratio), 
                                group=Group, color = Group)) + geom_line() +
      labs(title='Modification ratio (Relative to gene coverage)', subtitle=gene,
           x="Position", y ="Modification ratio (Relative to gene coverage)")+
      scale_fill_manual(values= c("#0000FF", "#FF0000"),
                        breaks=c(groups[1], groups[2]))+
      scale_x_discrete(limits = as.character(data_group1$Corrected_pos))+
      theme(axis.text=element_text(size=16, color="black"),
            axis.text.x = element_text(size=12,angle=90),
            axis.title= element_text(size=20), legend.text = element_text(size = 16),
            legend.title = element_text(size=20), plot.subtitle = element_text(size=19),
            plot.title =element_text(size=20))
    
    p
    file_name = paste0("Modification_ratio_by_pos_", gene, ".jpeg")
    ggsave(p, file=file_name, width = 35, 
           height = 25, units = "cm", path="../Results/Modification_ratio_plots/Comparison" )
    number_genes = number_genes+1
    ###FISHER TEST ####
    for(i in 1:length(positions)){
      if(positions[i] %in% data_group1$Corrected_pos){
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
      else{
        position_info[i] = 1
      }
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
adjusted_pvalues = p.adjust(aa_info$pvalue, method="BH", 
                            n=length(positions)*number_genes)
aa_info_final = data.frame(aa_info, adjusted_pvalues)
aa_info_final = aa_info_final[aa_info_final$adjusted_pvalues < 0.05,]


write.table(aa_info_final, file = "../Results/R_files/Modification_test.txt", 
            quote=FALSE, row.names = F)

#real_genes = substring(aa_info_final$gene, 1, nchar(aa_info_final$gene)-4)
unique_genes = unique(aa_info_final$gene)


for(aa in aas){
  id = grep(aa,unique_genes)
  cont = 1
  
  genes_aa = unique_genes[id]
  final_data = matrix(ncol=length(positions), nrow=length(genes_aa))
  rownames(final_data) = genes_aa
  colnames(final_data) = positions
  
  pvalue_list = final_data
  mod1_list = final_data
  mod2_list = final_data
  cov1_list = final_data
  cov2_list = final_data
  custom_text = final_data
  log2fc_list = final_data
  # Bucle for, in each isodecoder
  for(i in id){ 
    gene = unique_genes[i]
    data_group1 = read.table(file=paste0(dir,"/",groups[1],"_",gene, ".txt"),
                             header=T)
    data_group2 = read.table(file=paste0(dir,"/",groups[2],"_",gene, ".txt"), 
                             header=T)
  
    info_pvalue = c()
    info_mod1= c()
    info_mod2 = c()
    info_cov1 = c()
    info_cov2 = c()
    info_log2fc = c()

    for(pos in positions){
      nas = TRUE
      if(pos %in% aa_info_final$pos[aa_info_final$gene == gene]){
        cov1 = data_group1$Base_coverage[data_group1$Position==pos]
        cov2 = data_group2$Base_coverage[data_group2$Position==pos]
        if(cov1 > 50 | cov2 >50){
          mod1 = data_group1$Modification_ratio[data_group1$Position==pos] * 100
          mod2 = data_group2$Modification_ratio[data_group2$Position==pos] * 100
          if(mod1 >10 & mod2>10){
            nas = FALSE
  
            log2fc = (2*(mod2-mod1)/ (mod1+mod2))
            
            new_pvalue = (aa_info_final$adjusted_pvalues[aa_info_final$pos==pos & aa_info_final$gene==gene])
            new_pvalue = new_pvalue[1]  # Sometimes it returns the value repeated.
            new_pvalue = formatC(new_pvalue,format="e")
            info_pvalue = c(info_pvalue, new_pvalue)
            info_mod1 = c(info_mod1, mod1)
            info_mod2 = c(info_mod2, mod2)
            info_cov1 = c(info_cov1, cov1)
            info_cov2 = c(info_cov2, cov2)
            info_log2fc = c(info_log2fc, log2fc)
          }
          }
      }
      if(nas){
        info_pvalue = c(info_pvalue, NA)
        info_mod1 = c(info_mod1, NA)
        info_mod2 = c(info_mod2, NA)
        info_cov1 = c(info_cov1, NA)
        info_cov2 = c(info_cov2, NA)
        info_log2fc = c(info_log2fc, NA)
      }
      }
  
    
    final_data[cont,1:ncol(final_data)] = info_log2fc
    mod1_list[cont,1:ncol(final_data)] = info_mod1
    mod2_list[cont,1:ncol(final_data)] = info_mod2
    pvalue_list[cont,1:ncol(final_data)] = info_pvalue
    cov1_list[cont,1:ncol(final_data)] = info_cov1
    cov2_list[cont,1:ncol(final_data)] = info_cov2 
    log2fc_list[cont,1:ncol(final_data)] = info_log2fc
    cont = cont +1
  }
  if(nrow(final_data)>0){
    if(table(is.na(final_data)) != length(final_data)){
      setwd("../Results/Heatmaps")
      heatmap_file = paste0("Comparison_", aa, ".html")
      
      custom_text[] = paste0("Gene: ", rownames(final_data), "\n",
                             "Modification ", groups[1], " (%): ", mod1_list, "\n",
                             "Coverage ", groups[1], ": ", cov1_list, "\n",
                             "Modification ", groups[2], " (%): ", mod2_list, "\n",
                             "Coverage ", groups[2], ": ", cov2_list, "\n",
                             "(",groups[2], " - ", groups[1],")", "/mean: ", log2fc_list, "\n",
                             "Adjusted pvalue: ", pvalue_list)
      
      heatmap = heatmaply(final_data,  
                          colors= colorRampPalette(rev(brewer.pal(9, "RdBu"))),
                          plot_method = "plotly", 
                          limits=c(min(final_data, na.rm=TRUE),
                                   max(final_data,na.rm=TRUE)),
                          custom_hovertext=custom_text, Rowv = FALSE, Colv=FALSE, 
                          xlab="Position", ylab="Gene", column_text_angle=0, 
                          dendogram=FALSE, show_dendogram=c("FALSE", "FALSE"),
                          file=heatmap_file)
      setwd(dir_scripts)
    }
  } 
}

#out_table_x <- xtable(data_iso)
#print(out_table_x, type='html', file="../Results/DEG/Results_isotRNACP.html")

