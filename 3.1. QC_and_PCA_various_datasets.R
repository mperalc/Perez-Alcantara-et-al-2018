# QC and PCA RNA-seq data from various datasets

# This script creates Figure 1


currentDate <- Sys.Date() # to save date in name of output files

library(limma)
library(ggbiplot)
library(ggplot2)
library(reshape2) 
library(readr) 
library(org.Hs.eg.db)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
source("functions.R")
# functions


# filter counts function - compact version of what was done in script 1
# this script has an error: it didn't include the filter per stage and donor step that was done in QC 1. I checked and it would remove a further ~400 genes.
# shouldn't affect the PCA because the counts are very low and therefore don't drive the PCs.

filterCounts = function(counts){
  #   
  #   #get gene information for the Ensembl IDs
  # use archived mart for Grch38 (v88, march 2017) 
  
  ensembl38 = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "mar2017.archive.ensembl.org", path = "/biomart/martservice" ,dataset = "hsapiens_gene_ensembl")
  all_ensembl_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name'),
                            filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl38)
  rownames(all_ensembl_info) <- all_ensembl_info$ensembl_gene_id
  
  #check the geneID exists in Ensembl
  counts <- counts[all_ensembl_info$ensembl_gene_id,]
  
  #load the data into edgeR
  dge <- DGEList(counts = counts,genes = all_ensembl_info)
  
  #remove genes not expressed in at least one stage with CPM > 1
  isexpr <- rowSums(cpm(dge) > 1) >= 3   
  
  dge <- dge[isexpr, , keep.lib.sizes = FALSE]  # I had it previously (and always) on FALSE
  
  #remove genes that are not annotated as lincRNA or protein_coding
  dge <- dge[which(dge$genes$gene_biotype == "protein_coding" | dge$genes$gene_biotype == "lincRNA"),]
  #use only autosomal genes
  chr <- c(1:22)
  dge <- dge[which(dge$genes$chromosome_name %in% chr),]
  return(dge)
}

################## elements to be plotted

donors = c("Ad2.1","Ad3.1","Neo1.1")  # data comes from three donors
stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP","ENstage6", "ENstage7") # 8 stages
old_stages = c("iPSC","DE","PGT","PFG","PE","ENstage6") # 6 stages
old_donors = c("Ad2.1","Ad3.4")
Xie_stages = c("ES","DE","PGT","PFG","PE","late_PE","polyhormonal","matured_in_vivo")
Xie_donors = c(rep("Xie1",10),rep("Xie2",5),rep("Xie3",6),rep("Xie4",4))


##########
# load new differentiated cells' data

counts_new = read_tsv("counts.tsv") # load unfiltered conservative counts
counts_new = as.data.frame(counts_new)
rownames(counts_new) = counts_new[,1] # make gene id rownames
counts_new = counts_new[,c(-1:-2)] # take out first two columns
rownames(counts_new) <- gsub("(ENSG[0-9]+).*", "\\1", rownames(counts_new)) # cut everything after the dot
counts_new = orderColsByStages(counts_new) # order 


colnames(counts_new) = paste( rep(stage,each = 3),rep(donors,8),sep = "_" )  # rename columns

# old iPSC data Islets 2016
#75 bp read
counts_old = read.table("counts_old_trim75bp.tsv",header = T,check.names = F,row.names = 1)

counts_old <- counts_old[,13:ncol(counts_old)]
#cut the version number of the Ensembl IDs
rownames(counts_old) <- gsub("(ENSG[0-9]+).*", "\\1", rownames(counts_old))

colnames(counts_old) = paste(rep(old_stages,2),rep(old_donors,each = 6),sep = "_")  # change column names

counts_old = orderColsByStages(counts_old,stage = old_stages)


# Xie 2013

Xiecounts <- read.table("/Users/marta/Documents/WTCHG/DPhil/Data/Reference/Xie_2013_and_Nica/201215_stembancc.xie.nica.gene.counts",
                        header = T, check.names = F,row.names = 1)

#cut the version number of the Ensembl IDs
rownames(Xiecounts) <- gsub("(ENSG[0-9]+).*", "\\1", rownames(Xiecounts))
#get gene information for the Ensembl IDs
ensembl38 = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "mar2017.archive.ensembl.org", path = "/biomart/martservice" ,dataset = "hsapiens_gene_ensembl")
all_ensembl_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name'),
                          filters = 'ensembl_gene_id', values = rownames(Xiecounts), mart = ensembl38)
rownames(all_ensembl_info) <- all_ensembl_info$ensembl_gene_id

#check the geneID still exists in Ensembl
Xiecounts <- Xiecounts[all_ensembl_info$ensembl_gene_id,]

# remove unwanted samples:

Xiecounts = Xiecounts[,c(1,2,26:73)]


Xie <- c("DP1_ES_1","DP1_DE_1","DP1_GT_1","DP1_PF_1","DP1_PE_1",
         "DP1_ES_2","DP1_DE_2","DP1_GT_2","DP1_PF_2","DP1_PE_2",
         "DP3_ES","DP3_DE","DP3_GT","DP3_PF","DP3_PE",
         "DP3-4-5_late_PE_1","DP3-4-5_late_PE_2","DP3-4-5_late_PE_3",
         "DP3-4-5_poly_1","DP3-4-5_poly_2","DP3-4-5_poly_3", # PH
         "E2147_in-vivo_matured_1","E2147_in-vivo_matured_2","E2182_in-vivo_matured_1","E2182_in-vivo_matured_2") # FE
Xiecounts = Xiecounts[,Xie]  # sort

samples <- c(rep(donors,8),rep(old_donors,6),rep(Xie_donors))   # create the design matrix
samples = as.factor(samples)
stages <- c(rep(stage,each = 3),rep(old_stages,each = 2),
            rep(Xie_stages[1:5],3),rep(Xie_stages[6],3),rep(Xie_stages[7],3),rep(Xie_stages[8],4))
stages <- as.factor(stages)
study <- c(rep("new",24),rep("old",12),rep("Xie",25))
design <- model.matrix(~stages + samples + study)   # I'm not sure if I should group later stages than DE together.

combined_commongenes <- merge(counts_new, counts_old, by = 0, all = TRUE)  # merge by row names (by=0 or by="row.names")
rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
combined_commongenes <- combined_commongenes[,c(2:ncol(combined_commongenes))] # remove first column
combined_commongenes <- merge(combined_commongenes, Xiecounts, by = 0, all = TRUE)  # can only merge 2 dataframes at any time

combined_commongenes = na.omit(combined_commongenes)   # remove rows that contain NA values (reduces table to size of smallest table)
rownames(combined_commongenes) = combined_commongenes[,1] # make gene id rownames
combined_commongenes = combined_commongenes[,-1] # take out first column

filtered_combined_commongenes = filterCounts(combined_commongenes)   

filtered_combined_commongenes <- calcNormFactors(filtered_combined_commongenes)    # Calculate normalization factors. TMM by default

# save as an object for later analyses
#save(filtered_combined_commongenes, file = "/dge_old_new_Xie.xz" , compress="xz")

dim(filtered_combined_commongenes)


v <- voom(filtered_combined_commongenes,design,plot = F) # voom normalize the read counts

batch_corrected = removeBatchEffect(v$E,study)

######## remove unwanted experiments & samples

plotPcaDropExp = function(x, s = samples, st = stages){ # color: stage # shape: experiment -- drop samples
  pca1 <- prcomp(t(x), retx = TRUE)
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,Samples = samples,Stages = stages,Experiment = c(rep("Current",24),rep("Previous",12),rep("Xie",25)))
  levels(pcs$Stages) <- c(levels(pcs$Stages), c("EN","BLC","GT","PF","Late PE","Matured in vivo","Polyhormonal"))  # change name of stages, first increase levels of factor
  pcs[which(pcs$Stages == "ENstage6"),"Stages"] <- "EN"
  pcs[which(pcs$Stages == "ENstage7"),"Stages"] <- "BLC"
  pcs[which(pcs$Stages == "PGT"),"Stages"] <- "GT"
  pcs[which(pcs$Stages == "PFG"),"Stages"] <- "PF"
  pcs[which(pcs$Stages == "late_PE"),"Stages"] <- "Late PE"
  pcs[which(pcs$Stages == "matured_in_vivo"),"Stages"] <- "Matured in vivo"
  pcs[which(pcs$Stages == "polyhormonal"),"Stages"] <- "Polyhormonal"
  
  
  pcs <- droplevels(pcs)  # drop unused levels (old names)
  
  ################# remove unwanted experiments & samples after calculating the pcs
  
  pcs = pcs[which(!(pcs$Experiment == "Previous" & pcs$Stages %in% c("iPSC","DE","GT","PF","PE"))),]
  pcs = pcs[which(!(pcs$Experiment == "Xie" & pcs$Stages %in% c("ES","DE","GT","PF","PE","Late PE","Polyhormonal"))),]
  
  #################
  
  pcs$Stages <- ordered(pcs$Stages, levels = c("iPSC", "ES","DE", "GT", "PF", "PE", "Late PE","EP",
                                               "EN","Polyhormonal", "BLC","Matured in vivo") ) # order levels, for colours
  
  
  diaPalette <- c("#000000","#CADAE8", "#7883BA", "#755A91", "#CC85B1",  "#F4B8B0",
                  "#96665A", "#96165A","#6DA567")  # Diabetologia palette
  
  p <- ggplot(pcs,aes(x = PC1, y = PC2, color = Stages, shape = Experiment)) +
    geom_point(size = 5, aes(fill = Stages),
               #alpha = as.character(Experiment)),
               stroke = 1) +
    geom_point(size = 2.5,stroke = 1.5) +
    scale_fill_manual(values = c("Optimised" = diaPalette, "Previous" = diaPalette,"Xie" = "#FFFFFF")) +
    scale_colour_manual(values = diaPalette) +
    scale_shape_manual(values = c(16,17,22)) +
    xlab(paste0( "PC1: " ,percentVar[ 1 ],"% variance")) +
    ylab(paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  
  
  p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_blank(),
                    legend.key = element_blank(),
                    axis.title.y = element_text(  angle = 90, size = 24, vjust = 0.2),
                    axis.title.x = element_text(  size = 24, vjust = 0),
                    axis.text.x = element_text(  colour = "black", size = 24, vjust = 0.2, hjust = 1 ),
                    axis.text.y = element_text(  colour = "black",size = 24),
                    axis.ticks = element_line(colour = "black",size = 0.75),
                    axis.line = element_line(colour = "black",size = 0.30),
                    legend.text = element_text(  colour = "black",size = 24),
                    legend.title = element_text(  colour = "black",size = 24))
  
  plot(p)
  return(p)
}

tiff("new_diff_vs_selected_old_75bp_Xie_commongenes_PC1and2_batchcorrected_forDiabetologia_reduced.tiff", type = "cairo", antialias = "default",
     width = 10,height = 8,units = "in",res = 600,pointsize = 24,compression = "lzw")
 
p5_corrected_reduced = plotPcaDropExp(batch_corrected) 
 
dev.off()

