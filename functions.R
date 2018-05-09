# custom functions used in the analyses

orderColsByStages = function(counts,stage. = stage){
  nc = counts[ , grepl("iPSC", names(counts))]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in stage.[-1])  {         
    
    i = counts[ , grepl(s , names(counts))]
    nc = cbind(nc,i)
  }
  
  return(nc)
}


orderColsByDonors = function(x,donors. = donors){
  nc = x[ , grepl("Ad2.1" , colnames(x))]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in donors.[-1])  {         
    
    i = x[ , grepl(s , colnames(x))]
    nc = cbind(nc,i)
  }
  
  return(nc)
}


orderRowsByDonors = function(x,donors. = donors){
  nc = x[grep("Ad2.1", rownames(x)), ]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in donors.[-1])  {         
    
    i = x[grep(s, rownames(x)), ]
    nc = rbind(nc,i)
  }
  
  return(nc)
}

#function wrap so that plotMDS doesn't produce a plot when called:
plotMDS.invisible <- function(...){
  ff <- tempfile()
  png(filename = ff)
  res <- plotMDS(...)
  dev.off()
  unlink(ff)
  res
}


pretty_mds = function(v2){
  mds_p = plotMDS.invisible(v2$E,gene.selection = "pairwise")    # pairwise method (default)
  mds_c = plotMDS.invisible(v2$E,gene.selection = "common")      #common method 
  
  stages <- rep(stage,3) 
  # Rearrange data for ggplot
  
  # method: pairwise
  m_p = as.data.frame(mds_p$cmdscale.out)
  m_p <- cbind(m_p,sample = samples,stage = stages)
  colnames(m_p) = c("Dimension 1", "Dimension 2", "Samples","Stages")
  
  # method: common
  m_c = as.data.frame(mds_c$cmdscale.out)
  m_c <- cbind(m_c,sample = samples,stage = stages)
  colnames(m_c) = c("Dimension 1", "Dimension 2", "Samples","Stages")
  
  # plot pairwise
  
  mp = ggplot(m_p) +
    geom_point(aes(`Dimension 1` ,`Dimension 2`), size = 2, color = 'grey') +
    geom_label_repel(
      aes(`Dimension 1` , `Dimension 2`, fill = factor(Samples), label = Stages),
      fontface = 'bold', color = 'white',
      box.padding = unit(0.25, "lines"),
      point.padding = unit(0.25, "lines")
    ) +
    coord_fixed(ratio = 1.2) +
    
    theme_bw(base_size = 20) + 
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank() ) +
    theme(panel.border = element_rect()) +
    
    theme(legend.position = "bottom")  +
    ggtitle("MDS plot:pairwise method") +
    labs(fill = "Sample")
  
  
  # plot common
  
  mc = ggplot(m_c) +
    geom_point(aes(`Dimension 1` ,`Dimension 2`), size = 2, color = 'grey') +
    geom_label_repel(
      aes(`Dimension 1` , `Dimension 2`, fill = factor(Samples), label = Stages),
      fontface = 'bold', color = 'white',
      box.padding = unit(0.25, "lines"),
      point.padding = unit(0.25, "lines")
    ) +
    coord_fixed(ratio = 2) +  #fix x-y ratio
    
    theme_bw(base_size = 20) + 
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank() ) +
    theme(panel.border = element_rect()) +
    
    theme(legend.position = "bottom")  +
    ggtitle("MDS plot: common method") +
    labs(fill = "Sample")
  
  return(list(mp,mc))   #return both plots as a list
}



plot_sdc = function(x){
  sampleDists <- dist(t(x))   
  #This function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
  # by default, "euclidean"
  sampleDistMatrix <- as.matrix(sampleDists)
  stages <- rep(stage,3) 
  rownames(sampleDistMatrix) <- paste(stages, samples, sep = "-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  sdc = pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
  
  return(sdc)
}

plot_pca = function(x,s = Samples, st = Stages){
  pca1 <- prcomp(t(x), retx = TRUE)
  
  plot(pca1, type = "l") #variance vs first 10 components
  summary(pca1)          #importance of each component (important line is "proportion of variance")
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,Samples = Samples, Stages = Stages)
  pcs$Stages <- ordered(pcs$Stages, levels = Stages)
  
  diaPalette <- c("#CADAE8", "#7883BA", "#755A91", "#CC85B1", "#C15858", "#F4B8B0", "#96665A", "#6DA567")  # Diabetologia palette
  
  p <- ggplot(pcs, aes(PC1, PC2, colour = Stages, shape = Samples)) + 
    geom_point(size = 4) + xlab(paste0("PC1: " ,percentVar[ 1 ],"% variance")) + 
    ylab(paste0( "PC2: ",percentVar[ 2 ],"% variance" )) +
    scale_colour_manual(values = diaPalette)
  p <- p + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black"), 
                  legend.key = element_blank(),# legend.position = c(0.5,0.5),
                  axis.title.y = element_text(face = "bold", angle = 90, size = 16, vjust = 0.2),
                  axis.title.x = element_text(face = "bold", size = 16, vjust = 0),
                  axis.text.x = element_text(face = "bold", colour = "black", angle = 90, size = 16, vjust = 0.2, hjust = 1 ),
                  axis.text.y = element_text(face = "bold", colour = "black",size = 16),
                  axis.ticks = element_line(colour = "black"),
                  axis.line = element_line(colour = "black"),
                  legend.text = element_text(size = 14,face = "bold"),
                  legend.title = element_text(size = 16,face = "bold"))
  
  return(p)
}



make_long_plots = function(y, .tpm = tpm, .samples = samples, n = 9, .stage_2 = stage_2, name){
  
  plot_long = y[1:n,1:2]   # keep id and names to plot longitudinal variation of these top n genes 
  plot_long = tpm[match(plot_long$ensembl_gene_id,tpm$GeneID),]
  plot_long = orderColsByStages(plot_long)
  
  # melt data for ggplot2

  long = melt(plot_long, measure.vars = c(3:26))
  
  # rename stages and samples
  
  long$variable = rep(stage_2,each = n*3)
  colnames(long)[which(names(long) == "variable")] <- "stage"
  long$sample = rep(samples,each = n)
  long$stage <- factor(long$stage, levels = stage_2)
  long$sample = as.factor(long$sample)
  
  p <- ggplot(data = long, aes(x = stage, y = value,group = sample )) + ylab("transcript per million")  +
    geom_line(aes(col = sample), size = 1) + facet_wrap(~GeneName,ncol = 3,scales = "free") + ggtitle(paste("Longitudinal plots for top DE genes (contrasts):",name)) +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12))
 
  return(p)
}
