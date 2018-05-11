#
  #  PLOT EXPRESSION OF KEY GENES (IN TPM)
# This script creates ESM Figure 2
#   It needs data from TPM normalised counts - one can get that from the raw counts with featureCounts

# Load necessary libraries
library(readr)  # to read big tables fast
library(ggplot2) # to plot
library(reshape2)  # to modify dataframes for ggplot2
source("functions.R")

# Data from input columns
donor = c("Ad2.1", "Ad3.1", "Neo1.1")  # original samples, here called donor

stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "ENstage6", "ENstage7")  # 8 differentiation stages

# alternative names to plot:
stage_2 = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC")  #shortening EN names EN7= BLC (beta-like cells)

# rna-seq data
tpm = as.data.frame(read_tsv(tpm.tsv))  # read tpm file for plotting longitudinal tpm data

## ## ## ## CHANGES EVERY TIME ## ## ## ##
genes = c("POU5F1", "SOX17", "FOXA2", "NKX6-1", "NKX2-2", "NEUROG3", "MAFA", "MAFB",
          "INS", "ABCC8", "G6PC2", "GLP1R", "KCNJ11", "GCK", "SLC2A2", "GIPR", "SLC30A8")
# vector of genes to plot
output_name = "ESM_Figure_2"
output_type = "multiple_tiff"  # can be: 'png', 'tiff', 'multiple_png', 'multiple_tiff', 'pdf'
# number of columns for multiple plots
columns = 3
# colour palette
diaPalette <- c("#C15858", "#6DA567", "#7883BA")  # Diabetologia palette


########################################################
# QC and re-shaping of data for plotting

QC_and_reshaping = function() {
  # genes=sort(genes) # sort in alphabetical order
  plot_long = tpm[match(genes, tpm$GeneName), ]  # extracts from tpm data frame the rows that contain the genes of interest
  plot_long = na.omit(plot_long)  # remove NAs, in case there was not a match for a gene
  
  diff = setdiff(genes, tpm$GeneName)  # which ones in the genes' list are not in the table (probably have a different name)?
  diff
 
  plot_long = cbind(plot_long[c(1:2)], orderColsByStages(plot_long))
  gene_number = nrow(plot_long)  # how many genes to plot
  
  # melt data for ggplot2
  
  long = melt(plot_long, measure.vars = c(3:ncol(plot_long)))
  head(long)
  
  # rename stages and samples
  samples <- c(rep(donor, 8))
  long$variable = rep(stage_2, each = 3 * gene_number)  # sample size times number of genes
  colnames(long)[which(names(long) == "variable")] <- "stage"
  long$Sample = rep(samples, each = gene_number)
  long$stage <- factor(long$stage, levels = stage_2)
  long$Sample = as.factor(long$Sample)
  long$GeneName = factor(long$GeneName, levels = genes)
  return(long)
  
}  


if (output_type == "multiple_png" | output_type == "multiple_tiff" ) {
    ################ 1st part as QC script
    long = QC_and_reshaping()
    ###################### plot ########################
      
    p <- ggplot(data = long, aes(x = stage, y = value, group = Sample)) + 
      ylab("Expression (TPM)") + 
      expand_limits(y = 0) + 
      geom_hline(yintercept = 0, linetype = "dashed", size = 1, col = "#DCDCDC") + 
      geom_line(aes(linetype = Sample, col = Sample), size = 1) + 
      geom_point(size = 3, aes(shape = Sample, col = Sample)) + 
      scale_color_manual(values = diaPalette) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(size = 2), 
        axis.text.x = element_text(size = 16, face = "bold", angle = 45, vjust = 0.55), axis.text.y = element_text(size = 16, face = "bold"), axis.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), 
        legend.text = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        strip.background = element_blank(), strip.text = element_text(size = 16, face = "bold.italic")) +  # strip controls subfigure titles
      
      facet_wrap(~GeneName, scales = "free", ncol = columns)
    
    if (output_type == "multiple_png") {
      
      png(paste(output_name, ".png", sep = ""), type = "cairo", antialias = "default", width = 7, height = 7, units = "in", res = 600, pointsize = 13)
      print(p)
      dev.off()
    }
    
    if (output_type == "multiple_tiff") {
      tiff(paste(output_name, ".tiff", sep = ""), type = "cairo", compression = "lzw", antialias = "default", width = 14, height = 18, units = "in", res = 600, pointsize = 13)
      print(p)
      dev.off()
    }
  
} else {
  if (output_type == "pdf") {
   
    
    long = QC_and_reshaping()  # reshape for plotting
    head(long)
    
    #create plots first
    plot_list = list()
    
    for (i in unique(long$GeneName)) {
      long2 = long[long$GeneName == i, ]  # subset each gene
      p <- ggplot(data = long2, aes(x = stage, y = value, group = Sample)) +
        ggtitle(unique(long2$GeneName)) +
        xlab("Differentiation stages") +
        ylab("Expression [TPM]") +
        expand_limits(y = 0) +
        geom_hline(yintercept = 0, linetype = "dashed", size = 1, col = "#DCDCDC") + 
        geom_line(aes(linetype = Sample, col = Sample), size = 1) + 
        scale_colour_manual(values = diaPalette) +  # diabetologia pallete
        geom_point(size = 3,aes(shape = Sample,col = Sample)) +
        #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
        theme_bw() +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
              panel.border = element_rect(size = 2),axis.text = element_text(size = 12,face = "bold"),
              axis.title = element_text(size = 14,face = "bold"),plot.title = element_text(size = 16,face = "bold.italic"),
              legend.text = element_text(size = 11,face = "bold"),legend.title = element_text(size = 13,face = "bold"))
      
      plot_list[[i]] = p
      print(i)
    }
    
    # create progress bar
    library(tcltk)
    total = length(unique(long$GeneName))
    pb <- tkProgressBar(title = "progress bar", min = 0,
                        max = total, width = 300)
    
    # I have to open pdf connection before looping so that it saves one gene on each page
    somePDFPath = paste(output_name,".pdf",sep = "")
    pdf(file = somePDFPath)
    j = 1
    for (i in unique(long$GeneName)) {
      print(plot_list[[i]])
      Sys.sleep(0.1)
      setTkProgressBar(pb, j, label = paste( round(j/total*100, 0),
                                           "% done"))
      j = j + 1
    }
    close(pb)
    dev.off()
    
  }else{
    
    for (genes in genes) { # individually gene by gene
      
      long = QC_and_reshaping()  # reshape for plotting
      head(long)
      
      p <- ggplot(data = long2, aes(x = stage, y = value, group = Sample)) +
        ggtitle(unique(long2$GeneName)) +
        xlab("Differentiation stages") +
        ylab("Expression [TPM]") +
        expand_limits(y = 0) +
        geom_hline(yintercept = 0, linetype = "dashed", size = 1, col = "#DCDCDC") + 
        geom_line(aes(linetype = Sample, col = Sample), size = 1) + 
        scale_colour_manual(values = diaPalette) +  # diabetologia pallete
        geom_point(size = 3,aes(shape = Sample,col = Sample)) +
        #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
        theme_bw() +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
              panel.border = element_rect(size = 2),axis.text = element_text(size = 12,face = "bold"),
              axis.title = element_text(size = 14,face = "bold"),plot.title = element_text(size = 16,face = "bold.italic"),
              legend.text = element_text(size = 11,face = "bold"),legend.title = element_text(size = 13,face = "bold"))
      
      if (output_type == "png") {
        png(paste(output_name,".png",sep = ""), type = "cairo",
            width = 8,height = 5,units = "in",res = 300,pointsize = 12)
        print(p)
        dev.off()
        
      }
      if (output_type == "tiff") {
        tiff(paste(output_name,".tiff",sep = ""), type = "cairo", compression = "lzw", antialias = "default",
             width = 8,height = 5,units = "in",res = 1000,pointsize = 13)
        print(p)
        dev.off()
      }
    }
  }
  
}  
