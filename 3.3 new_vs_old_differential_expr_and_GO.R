# compare old beta-cell differentiation model (pre-2016) to new one (2016-2017).
# Differential expression analysis of old vs new data (ESM Table 7) and gene ontology analysis

currentDate <- Sys.Date() # to save date in name of output files

library(limma)
library(ggbiplot)
library(ggplot2)
library(reshape2)  # to modify dataframes for ggplot2
library(readr) # fast reading of large files
library(org.Hs.eg.db)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(GOstats)
library(qvalue)
source("functions.R")

## create folder for results under the current directory
dir.create(file.path(getwd(), "new_vs_old_differential_expression_results"))


# filter counts

# filter counts function - compact version of what was done in script 1
# this script has an error: it didn't include the filter per stage and donor step that was done in QC 1. 

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


# plot PCA
plotPCA = function(x,s = samples, st = stages) {
  pca1 <- prcomp(t(x), retx = TRUE)
  
  percentVar <- (pca1$sdev)^2/sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs, Samples = samples, Stages = stages, Experiment = c(rep("New", 24), rep("Old", 12)))
  levels(pcs$Stages) <- c(levels(pcs$Stages), c("EN", "BLC", "GT", "PF"))  # change name of stages, first increase levels of factor
  pcs[which(pcs$Stages == "ENstage6"), "Stages"] <- "EN"
  pcs[which(pcs$Stages == "ENstage7"), "Stages"] <- "BLC"
  pcs[which(pcs$Stages == "PGT"), "Stages"] <- "GT"
  pcs[which(pcs$Stages == "PFG"), "Stages"] <- "PF"
  
  pcs <- droplevels(pcs)  # drop unused levels (old names)
  
  pcs$Stages <- ordered(pcs$Stages, levels = c("iPSC", "DE", "GT", "PF", "PE", "EP","EN", "BLC") ) # order levels, for colours
  
  diaPalette <- c("#CADAE8", "#7883BA", "#755A91", "#CC85B1", "#C15858", "#F4B8B0", "#96665A", "#6DA567")  # Diabetologia palette
  
  p <- ggplot(pcs, aes(x = PC1, y = PC2, color = Stages, shape = Samples)) +
    geom_point(size = 4, aes(fill = Stages, alpha = as.character(Experiment)), stroke = 1) +
    geom_point(size = 2.5, stroke = 1.5) +      
    scale_shape_manual(values = c(22,24,25,21)) +
    scale_alpha_manual(values = c(Old = 0, New = diaPalette), name = "Experiment", guide = "none") + # guide="none" takes out legend for this alpha
    scale_colour_manual(values = diaPalette) +
    scale_fill_manual(values = diaPalette) +
    xlab(paste0( "PC1:" ,percentVar[ 1 ],"% variance")) + 
    ylab(paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  
  
  p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black"), 
                    legend.key = element_blank(),# legend.position = c(0.5,0.5),
                    axis.title.y = element_text(face = "bold", angle = 90, size = 16, vjust = 0.2),
                    axis.title.x = element_text(face = "bold", size = 16, vjust = 0),
                    axis.text.x = element_text(face = "bold", colour = "black", angle = 90, size = 16, vjust = 0.2, hjust = 1),
                    axis.text.y = element_text(face = "bold", colour = "black", size = 16),
                    axis.ticks = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    legend.text = element_text(face = "bold", colour = "black",size = 14),
                    legend.title = element_text(face = "bold", colour = "black",size = 16))
  
  return(p)
}

########## Data filtering ##############

# load new data
counts = read_tsv("counts.tsv") # load unfiltered conservative counts
counts = as.data.frame(counts)
rownames(counts) = counts[, 1] # make gene id rownames
counts = counts[, c(-1:-2)] # take out first two columns
rownames(counts) = gsub("(ENSG[0-9]+).*", "\\1", rownames(counts)) # cut everything after the dot

donors = c("Ad2.1", "Ad3.1", "Neo1.1")  # data comes from three donors
stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "ENstage6", "ENstage7") # 8 stages
old_stages = c("iPSC", "DE", "PGT", "PFG", "PE", "ENstage6") # 6 stages
old_donors = c("Ad2.1", "Ad3.4")

counts = orderColsByStages(counts) # order
colnames(counts) = paste(rep(stage, each = 3), rep(donors, 8), sep = "_")  # rename columns

# load old differentiation data
#75 bp read
old_diff = read.table(
  "counts_old_trim75bp.tsv",
  header = T,
  check.names = F,
  row.names = 1
)

old_diff <- old_diff[, 13:ncol(old_diff)]
#cut the version number of the Ensembl IDs
rownames(old_diff) = gsub("(ENSG[0-9]+).*", "\\1", rownames(old_diff))
colnames(old_diff) = paste(rep(old_stages, 2), rep(old_donors, each = 6), sep = "_")  # change column names
old_diff = orderColsByStages(old_diff, stage = old_stages)

samples <- c(rep(donors, 8), rep(old_donors, 6))
samples = as.factor(samples)
stages <- c(rep(stage, each = 3), rep(old_stages, each = 2))
stages <- as.factor(stages)
study <- c(rep("new", 24), rep("old", 12))
design <-  model.matrix( ~ stages + samples + study)   


combined_commongenes <-  merge(counts, old_diff, by = 0, all = TRUE)  # merge by row names 
combined_commongenes = na.omit(combined_commongenes)   # remove rows that contain NA values 
rownames(combined_commongenes) = combined_commongenes[, 1] # make gene id rownames
combined_commongenes = combined_commongenes[, -1] # take out first column


filtered_combined_commongenes = filterCounts(combined_commongenes)   
nrow(filtered_combined_commongenes)

filtered_combined_commongenes <-  calcNormFactors(filtered_combined_commongenes)    # Calculate normalization factors. TMM by default
# save as an object for later analyses
#save(filtered_combined_commongenes, file = "dge_old_and_new_filtered_75bp.xz" , compress="xz")

v <-  voom(filtered_combined_commongenes, design, plot = F) # voom normalize the read counts

# PCA
p = plotPCA(v$E)
p

######### differential expression analysis old vs new data ###########
# for each stage, throughout all shared stages (iPSC to EN6)

#load(file="dge_old_and_new_filtered_75bp.xz",verbose=TRUE)  #loading the dge object for conservative counts

design <- model.matrix( ~ stages + samples)
experimentnew <- c(rep(1, 24), rep(0, 12))  # base level: old experiment
design = cbind(design, experimentnew)  # add to design

# differential expression analysis old vs new, 
# reorder original design matrix, and take out latest stage to make looping through column names easier:

design = cbind(design[, c(1, 12, 5)], c(rep(0, 3), rep(1, 3), rep(0, 30)), design[, c(8, 7, 6, 4, 2, 9, 10, 11)])
colnames(design)[4] <- "stagesDE"

diff_exp_allstages <- list()
diff_exp_allstages_nominal <- list()
for (s in old_stages) {
  v2 = voom(filtered_combined_commongenes,
            design = design,
            save.plot = TRUE)
  rownames(design) = colnames(v2$E)   # to select rows more easily
  v3 = v2[, grepl(s , colnames(v2))]  # take only columns from selected stage
  
  design3 = design[grepl(s , rownames(design)), c(1, 2)]  # same with design matrix
  
  # filter per stage to make sure I'm not comparing leftover genes with 
  # keep genes that have cpm>1 from voom object in all samples in either experiment old or new in the stage we are checking
  cpm_stage = cpm(filtered_combined_commongenes)  # calculate cpm
  cpm_stage = cpm_stage[, grepl(s , colnames(cpm_stage))]  # select columns from stage
  
  cpm_stage[cpm_stage < 1] = 0  # change to binary to make the selection easier
  cpm_stage[cpm_stage >= 1] = 1
  keep_genes = names(which(rowSums(cpm_stage[, c(1:3)]) == 3 |
                             rowSums(cpm_stage[, c(4:5)]) == 2))  # select genes that pass filter
  v3 = v3[keep_genes, ]
  
  ###
  fit3 <- lmFit(v3, design3)
  fit3 <- eBayes(fit3)
  
  diff_exp3 = topTable(
    fit3,
    coef = ncol(design3),
    sort.by = "none",
    number = nrow(fit3$coefficients)
  )
  diff_exp_allstages_nominal[[s]] = diff_exp3[which(diff_exp3$P.Value <
                                                      0.01 & (diff_exp3$logFC > 1 | diff_exp3$logFC < (-1))), ]
  diff_exp_allstages[[s]] = diff_exp3[which(diff_exp3$adj.P.Val < 0.01 &
                                              (diff_exp3$logFC > 1 |
                                                 diff_exp3$logFC < (-1))), ]  #pval<0.01 and logFC>1
  write.csv(
    diff_exp_allstages[[s]],
    file = paste(
      "./new_vs_old_differential_expression_results/",
      currentDate,
      "_diff_exp_new_vs_old_protocol_",
      s,
      "_stage.csv",
      sep = ""
    ),
    row.names = F
  )
  #  nrow-ing each table from the list shows a slight mistake in the paper: in the last stage there are 2094 diff expressed genes, 1 less than 
  # stated in the paper
}

############### GO analysis ###############################

stage = c("iPSC", "DE", "PGT", "PFG", "PE", "ENstage6")

# background list
# load(
#   "dge_old_and_new_filtered_75bp.xz"
# )  

# test lists
sig_stages = list()

for (s in stage) {
  # read in DEA data for each stage
  sig_stages[[s]] <-
    read.csv(
      file =  paste(
        "./new_vs_old_differential_expression_results/",
        currentDate,
        "_diff_exp_new_vs_old_protocol_",
        s,
        "_stage.csv",
        sep = ""
      ),
      header = T
    )
  
  sig_stages[[s]] = sig_stages[[s]][order(sig_stages[[s]]$adj.P.Val, decreasing =
                                            F), ] # order by adjusted p.value
}

##get the entrez IDs for the genes of interest. Takes Ensembl gene ids from my data


ensembl38 = useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  host = "mar2017.archive.ensembl.org",
  path = "/biomart/martservice" ,
  dataset = "hsapiens_gene_ensembl"
)
all_ensembl_to_entrez = getBM(
    attributes = c('ensembl_gene_id', 'entrezgene'),
    filters = 'ensembl_gene_id',
    values = rownames(filtered_combined_commongenes$genes),
    mart = ensembl38
  )  # all genes from dge object


entrez_object <- org.Hs.egGO
mapped_genes <- mappedkeys(entrez_object) # all entrez keys that have GO id
all_with_go <- unique(intersect(all_ensembl_to_entrez$entrezgene,mapped_genes)) # all entrez ids that have GO id from dge object


## for genes with higher diff expr in old experiment:

sig_stages_neg <- list()
for (s in stage) {
  sig_stages_neg[[s]] <- sig_stages[[s]][which(sig_stages[[s]]$logFC < 1), ]
}

sig_genes <- list()

for (i in stage) {
  sig_genes[[i]] = getBM(
      attributes = c('ensembl_gene_id', 'entrezgene'),
      filters = 'ensembl_gene_id',
      values = sig_stages_neg[[i]]$ensembl_gene_id,
      mart = ensembl
    )
}

sig_genes_with_go = list()
for (i in stage) {
  sig_genes_with_go[[i]] = unique(intersect(sig_genes[[i]]$entrezgene, mapped_genes))
}

#perform GO enrichment using GOstats


GO_list = list()

for (i in stage) {
  params = new(
    'GOHyperGParams',
    geneIds = sig_genes_with_go[[i]],
    # list of genes I'm testing
    universeGeneIds = all_with_go,
    # list of genes I used in my DE analysis (all my genes with RNA-seq data)
    ontology = 'BP',
    # biological process
    pvalueCutoff = 1,
    # p value=1 to be able to calculate FDR adjusted p values later (q-values)
    conditional = T,
    testDirection = 'over',
    annotation = "org.Hs.eg.db"
  )
  hgOver = hyperGTest(params)
  hgOver
  
  result = summary(hgOver)
  result <- result[which(result$Size > 1 & result$Count > 1),] 
  
  
  # # ## calculate adjusted pvalues (BH)
  result <- cbind(result, adjust_pval = p.adjust(result$Pvalue, "BH"))
  result <- result[which(result$adjust_pval < 0.05), ] # saving GO terms with significant adjusted p-values
  
  GO_list[[i]] <- result
  
}

for (i in stage) {
  # change name according to method used to calculate adjusted p-vals and logFC>1 or 0
  write.csv(
    file = paste(
      "./new_vs_old_differential_expression_results/",
      currentDate,
      "_GOstats_GO_BP_",
      i,
      "expr_higherinold.csv",
      sep = ""
    ),
    GO_list[[i]],
    row.names = F
  )
}

# higher in new:

sig_stages_pos <- list()
for (s in stage) {
  sig_stages_pos[[s]] <-
    sig_stages[[s]][which(sig_stages[[s]]$logFC > 1), ]
}

sig_genes <- list()
for (i in stage) {
  sig_genes[[i]] <-
    getBM(
      attributes = c('ensembl_gene_id', 'entrezgene'),
      filters = 'ensembl_gene_id',
      values = sig_stages_pos[[i]]$ensembl_gene_id,
      mart = ensembl38
    )
}

sig_genes_with_go <- list()
for (i in stage) {
  sig_genes_with_go[[i]] = unique(intersect(sig_genes[[i]]$entrezgene, mapped_genes))
}

GO_list <- list()
for (i in stage) {
  params <- new(
    'GOHyperGParams',
    geneIds = sig_genes_with_go[[i]],
    # list of genes I'm testing
    universeGeneIds = all_with_go,
    # list of genes I used in my DE analysis (all my genes with RNA-seq data)
    ontology = 'BP',
    # biological process
    pvalueCutoff = 1,
    # p value=1 to be able to calculate FDR adjusted p values later (q-values)
    conditional = T,
    testDirection = 'over',
    annotation = "org.Hs.eg.db"
  )
  hgOver <- hyperGTest(params)
  hgOver
  
  result <- summary(hgOver)
  
  
  result <- result[which(result$Size > 1 & result$Count > 1), ]
  
  # # ## calculate adjusted pvalues (BH)
  result <- cbind(result, adjust_pval = p.adjust(result$Pvalue, "BH"))
  result <-
    result[which(result$adjust_pval < 0.05), ] # saving GO terms with significant adjusted p-values
  
  GO_list[[i]] <- result
  
}


for (i in stage) {
  # change name according to method used to calculate adjusted p-vals and logFC>1 or 0
  write.csv(
    file = paste(
      "./new_vs_old_differential_expression_results/",
      currentDate,
      "_GOstats_GO_BP_",
      i,
      "expr_higherinnew.csv",
      sep = ""
    ),
    GO_list[[i]],
    row.names = F
  )
}
