#! /usr/bin/Rscript

list.of.packages <- c("GenomicAlignments","GenomicFeatures", "Rsamtools","edgeR","DESeq2",
                      "vsn", "magick", "Glimma","pheatmap", "tidyverse", "readxl",
                      "data.table","plyr") 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  install.packages(new.packages, repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages({
library(Rsubread, quietly = TRUE, warn.conflicts = FALSE)
library( "GenomicAlignments", quietly = TRUE, warn.conflicts = FALSE )
library( "GenomicFeatures", quietly = TRUE, warn.conflicts = FALSE )
library( "Rsamtools", quietly = TRUE, warn.conflicts = FALSE )
library("edgeR", quietly = TRUE, warn.conflicts = FALSE)
library("DESeq2", quietly = TRUE, warn.conflicts = FALSE)
library("vsn", quietly = TRUE, warn.conflicts = FALSE)
library("magick", quietly = TRUE, warn.conflicts = FALSE)
library("Glimma", quietly = TRUE, warn.conflicts = FALSE)
library("pheatmap", quietly = TRUE, warn.conflicts = FALSE)
library("tidyverse")
library("readxl")
library("plyr")
library("data.table")
})

sample_data <- read.table("../Fastq_downloaded/sample_data.txt",
                               header=T, stringsAsFactors = T)

#Read proportions table
counts <- read.table("../Results/R_files/Counts/Counts_by_gene_total_for_DESeq2.txt", 
                          header=TRUE, quote="\"", row.names=1)


counts = data.frame(counts, TRNA=rownames(counts))

#obtain prop
counts <- data.frame(counts,do.call(rbind,strsplit(as.character(counts$TRNA),"-")))
counts[ ,c('V4', 'X1','X4','X5')] <- list(NULL)
colnames(counts)[length(counts)]<- "Cod"
colnames(counts)[length(counts)-1]<- "Aa"

counts = counts[counts$Cod!="*",]
counts = counts[order(counts$Cod),]
cod_list <- unique(counts$Cod)

all_prop <- data.frame(matrix(ncol = ncol(counts), nrow = 0))
colnames(all_prop)<- colnames(counts)

for (cod in cod_list) {
  trna_cod <- counts[counts$Cod == cod, ]
  sum <- apply(trna_cod[1:(ncol(trna_cod)-3)], 2, function(x) sum(x))
  trna_cod = as.data.frame(t(apply(trna_cod[1:(ncol(trna_cod)-3)], 1, 
                                         function(x) x/sum)))
  all_prop <- rbind.fill(all_prop,trna_cod)
}

rownames(all_prop) = rownames(counts)
all_prop$TRNA = counts$TRNA
all_prop$Cod = counts$Cod
all_prop$Aa = counts$Aa


proportions = all_prop[c(1:(ncol(all_prop)-3))]



options(scipen=999)


conditions1 = sample_data$Condition == levels(sample_data$Condition)[1]
conditions2 = sample_data$Condition == levels(sample_data$Condition)[2]
# Apply Wilcoxon test (a nonparametric statistical test that compares two paired groups)
is0 = apply(proportions, 1, function(x) sum(x)==0)
proportions = na.omit(proportions[!is0,])
test <- apply(proportions,1, function(a) 
  wilcox.test(x=a[conditions1],y=a[conditions2], exact=FALSE)$p.value)
test = as.data.frame(test)
test_adjust = as.data.frame(p.adjust(test, method="BH"))

write.table(test, file = "../Results/R_files/p-val_iso-trna-cp.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

options(scipen = 999)

# No vale Glimma.

group1 = which(sample_data$ID== sample_data$ID[
  sample_data$Condition ==  levels(sample_data$Condition)[1] ] )
group2 = which(sample_data$ID== sample_data$ID[
  sample_data$Condition ==  levels(sample_data$Condition)[2] ] )

unique_aa = unique(all_prop$Aa)

if (!file.exists("../Results/Plots/isotRNACP")){
  dir.create("../Results/Plots/isotRNACP", recursive=TRUE)
}
for(aa in unique_aa){
  data = all_prop[all_prop$Aa == aa,]
  first_group = apply(data[group1],1, function(x) mean(x))
  second_group = apply(data[group2],1, function(x) mean(x))
  means = as.numeric(c(first_group, second_group))
  group = c(rep(levels(sample_data$Condition)[1], length(first_group)),
            rep(levels(sample_data$Condition)[2], length(second_group)))
  
  data_final = data.frame(cbind(means, group, tRNA =rep(data$TRNA,2), 
                                Aa=rep(data$Aa,2), Cod=rep(data$Cod,2)))
  data_final$means = as.numeric(data_final$means)
  p= ggplot(data = data_final, aes(x=tRNA, y=means, fill=group))+ 
    scale_fill_manual(values=c("dodgerblue","olivedrab3"))+
    geom_bar(stat="identity", position=position_dodge2()) + 
    labs(x="tRNA family",y = "Proportion", fill="Group")+ 
    theme(text = element_text(size=10), 
          axis.text.x = element_text(angle=90, color="black", vjust=0.5), 
          axis.text.y = element_text(color="black"), 
          panel.background =element_rect(fill = "snow2", colour = "snow2",
                                    size = 0.5, linetype = "solid"),
          strip.text.x = element_text(face="bold")) + 
    facet_grid(~Cod, scales="free", space="free")
  file_name = paste0("../Results/Plots/isotRNACP/Plots_iso_", aa, ".jpeg")
  ggsave(plot=p, filename=file_name,
         width = 20, height = 10, units = "cm")
  
}


####MEAN


mean_group1 = apply(proportions,1, function(x) mean(x[group1]))
mean_group2 = apply(proportions,1, function(x) mean(x[group2]))
diff_expr = mean_group1 - mean_group2

mean_data = data.frame(mean_group1, mean_group2)
colnames(mean_data) = c(levels(sampleInfo_total$Condition))

jpeg("../Results/Plots/Heatmap_mean_proportions.jpeg")
p = pheatmap(mean_data, cluster_rows=T, show_rownames=T, show_colnames = T,
         cluster_cols=F,fontsize_row=0.5, clustering_method="ward.D")
res_data_mean <- as.data.frame(mean_data)
write.table(res_data_mean, file = "../Results/R_files/Counts_mean_heatmap.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
dev.off()



###
mean_data <- transform(proportions, Control=(proportions$RC01_Control+proportions$RC02_Control+proportions$RC03_Control)/3)

mean_data <- transform(mean_data, Treated=(proportions$RC05_Treated+proportions$RC06_Treated+proportions$RC08_Treated)/3)
mean_data <- subset(mean_data, select=c(Control,Treated))

pdf("../Results/heatmap_iso_all_CRG_mean.pdf")
pheatmap(mean_data, cluster_rows=T, show_rownames=T, show_colnames = T,
         cluster_cols=FALSE,fontsize_row=0.5, clustering_method="ward.D")

#Differentially expressed genes (TRUE = p-adj < 0.05 | FALSE = p-adj > 0.05)

test_res05 <- test < 0.05
sum <- summary(test_res05)
sum
names(sum)[names(sum) == "FALSE"] <- "[p-adj > 0.05]"
names(sum)[names(sum) == "TRUE"] <- "[p-adj < 0.05]"
names(sum)[names(sum) == "Mode"] <- ""
names(sum)[names(sum) == "logical"] <- "Number of genes"

