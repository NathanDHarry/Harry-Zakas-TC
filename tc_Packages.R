####
####  Packages for Blueberry RNAseq analysis
####

# detach any currently loaded packages
lapply(names(sessionInfo()$otherPkgs), function(pkgs)
  detach(paste0('package:', pkgs), character.only = T, unload = T, force = T))

##
## Package Installation ----
##

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
# Package installation function for packages from CRAN
using.cran<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}
# Package installation function for packages from bioconda
using.biocon<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    BiocManager::install(need)
    lapply(need,require,character.only=TRUE)
  }
}

### CRAN required packages ----
req.package.cran <- c("pheatmap", "plyr", "tidyverse", "ggthemes", "gridExtra", "scales", "ape")

### Bioconda required packages ----
req.package.biocon <- c("DESeq2", "tximport", "GenomicFeatures", "Mfuzz", "topGO")

## Execute the functions to install or load packages
using.cran(req.package.cran)
using.biocon(req.package.biocon)

###
### Settings (variables) ----
###

code.Ambiguous <- c(1,2,10,12,20,21,100,101,120,121,200,202,210,212)
code.Additive <- c(112, 221)
code.Conserved <- c(0)
code.LDominant <- c(102, 201)
code.PDominant <- c(110, 220)
code.OverDominant <- c(22, 122, 222)
code.UnderDominant <- c(11, 111, 211)


##
## Custom Functions ----
##

plot.data_PCA <- function(object, intgroup = "Genotype", ntop = 500, returnData = FALSE){
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = T)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))){
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = F])
  intgroup.df <- lapply(intgroup.df, factor)
  d <- data.frame("PC1" = pca$x[, 1], "PC2" = pca$x[, 2], "PC3" = pca$x[, 3], "PC4" = pca$x[, 4], "PC5" = pca$x[, 5], intgroup.df, "name" = colData(object)$sampleID)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:4]
    return(d)
  }
}

plot.PCA <- function(PCAData, PCa, PCb, colo, shape = NULL){
  Pa <- as.character(paste("PC", PCa, sep=""))
  Pb <- as.character(paste("PC", PCb, sep=""))
  colo <- as.character(colo)
  colors.select <- list("Genotype" = scale_color_manual(values = colors.streb[(names(colors.streb) %in% PCAData$Genotype)]),
                        "seq_Run" = scale_color_manual(values = c("tc1" = "royalblue", "tc2" = "darkred")))
  
  both.min <- min(PCAData[,PCa], PCAData[,PCb]) - 5
  both.max <- max(PCAData[,PCa], PCAData[,PCb]) + 5
  
  ggplot(PCAData, aes_string(x = Pa, y = Pb, color = colo, label = "name", shape = shape)) +
    geom_point(size = 1, stroke = 1.5, alpha = 0.75) +
    xlab(paste0(paste(Pa, ": ", sep = ""), round(attr(PCAData, "percentVar")[as.numeric(PCa)] * 100), "% variance")) +
    ylab(paste0(paste(Pb, ": ", sep = ""), round(attr(PCAData, "percentVar")[as.numeric(PCb)] * 100), "% variance")) +
    colors.select[[colo]] +
    #scale_y_continuous(limits = c((min(PCAData[,PCb]) - 5), (max(PCAData[,PCb]) + 5))) +
    scale_y_continuous(limits = c(both.min, both.max)) +
    #scale_x_continuous(limits = c((min(PCAData[,PCa]) - 5), (max(PCAData[,PCa]) + 5))) +
    scale_x_continuous(limits = c(both.min, both.max)) +
    coord_fixed() +
    labs(title = "Gene Expression", color = "black") +
    theme(axis.text = element_text(size = 15), axis.title.x = element_text(margin = margin(t = 18,)),
          axis.title = element_text(size = 18), plot.title = element_text(size = 20), legend.text = element_text(size = 10), legend.title = element_text(size = 15)) +
    guides(color = guide_legend(gsub("\"","",deparse(substitute(colo)))), shape = guide_legend(gsub("\"","",deparse(substitute(shape)))))
}



# Function to plot individual genes over time in P and L samples
plot.geneTC <- function(gene_ID, ds.object=dds, metadata=data.meta.trim){
  time_points <- c(1, 2, 3, 4, 5, 6)
  names(time_points) <- levels(metadata$TissueType)
  
  tc_gene.df <- data.frame("norm_counts" = DESeq2::counts(ds.object, normalized = T)[gene_ID, metadata$sampleID[grepl("^LL$|^PP$", metadata$Genotype, perl=T)]])
  tc_gene.df$TissueType <- time_points[metadata$TissueType[grepl("^LL$|^PP$", metadata$Genotype, perl=T)]]
  tc_gene.df$Genotype <- metadata$Genotype[grepl("^LL$|^PP$", metadata$Genotype, perl=T)]
  
  d <- ggplot(data = tc_gene.df, mapping = aes(x = TissueType, y = norm_counts, col = Genotype)) +
    #geom_point() + 
    geom_smooth(method = "loess", formula = y ~ x, se = F, linewidth = 1.5) +
    scale_x_continuous(name = "Developmental Stage", 
                       breaks = time_points, 
                       labels = levels(metadata$TissueType),
                       limits = c(0.9, 6.1)) +
    scale_y_continuous(trans = 'log10') +
    labs(x = "Developmental Time", y = "Normalized Expression Counts (log10)", caption = paste0(GeneID, ": \n", annot)) +
    scale_color_manual(values = colors.streb[1:2])
  
  return(d)
}


# Function to plot individual genes over time in P and L samples with ridge-plots
ridgeplot.geneTC <- function(gene_ID, ds.object=dds, metadata=data.meta.trim){
  time_points <- c(1, 2, 3, 4, 5, 6)
  names(time_points) <- levels(metadata$TissueType)
  
  tc_gene.df <- data.frame("norm_counts" = DESeq2::counts(ds.object, normalized = T)[gene_ID, metadata$sampleID[grepl("^LL$|^PP$", metadata$Genotype, perl=T)]])
  tc_gene.df$TissueType <- time_points[metadata$TissueType[grepl("^LL$|^PP$", metadata$Genotype, perl=T)]]
  tc_gene.df$Genotype <- metadata$Genotype[grepl("^LL$|^PP$", metadata$Genotype, perl=T)]
  tc_gene.df$Genotype <- factor(tc_gene.df$Genotype, levels = c("PP", "LL"))
  dta <- aggregate(norm_counts ~ Genotype:TissueType, tc_gene.df, mean)
  
  
  d <- ggplot(data = dta, mapping = aes(x = `TissueType`, y = `Genotype`, height = `norm_counts`, fill = `Genotype`)) +
    geom_density_ridges(stat = "identity", scale = 1.25) +
    scale_fill_manual(values = colors.streb[1:2]) +
    scale_y_discrete(expand = c(0.05, 0)) +  
    scale_x_continuous(expand = c(0, 0), name = "Time Point", 
                       breaks = time_points, 
                       labels = levels(metadata$TissueType),
                       limits = c(0.9, 6.1)) +
    labs(y = "Average Gene Expression")
  
  return(d)
}

# DAFS qualitative thresholding function (requries mclust and earth packages)
DAFS_qcut <- function(data) {
  
  out <- apply(data,1,function(x) all(x==0))
  if(length(out[out=="TRUE"])>0) data <- data[-which(out=="TRUE"),]
  
  
  #set vector for cutoff values
  cutv <- rep(0,0)
  
  for (i in 1:ncol(data)) {
    
    #specify array and remove 0 counts
    xx <- data[,i]
    xx <- xx[-which(xx==0)]
    
    #take log2 of data
    log2xx <- log2(xx)
    dlog2 <- data.frame(LogC=log2xx)
    
    #vector to store Kolmogorov Smirnov distance statistics
    vv <- rep(0,0)
    
    #select start point
    start <- length(log2xx[log2xx==min(log2xx)])/length(log2xx)
    
    #set sequence
    s <- seq(round(start,2),0.5,by=0.005)
    
    #loop through cuts of the data to determine targeted K-S statistic
    for(q in s) {
      
      #select data greater than a quantile and run Mclust on that data to determine theoretical distribution
      d <- log2xx[which(log2xx>quantile(log2xx,q,na.rm=T))]
      out <- Mclust(d,G=1)
      ks <- ks.test(d,"pnorm",out$parameter$mean,
                    out$parameter$variance$sigmasq)
      vv <- c(vv,ks$statistic)
    }
    
    #determine first left-most local minima
    qq <- earth(s,vv,thresh=0.005)
    
    #save suggested cut
    cutv <- c(cutv,min(qq$cuts[qq$cuts>0]))
  }
  
  names(cutv) <- colnames(data)
  return(cutv)
}

## Clean up environment
rm(using.cran, using.biocon, req.package.cran, req.package.biocon)


# Plot MFuzz clusters
#png(filename = "CLUSTERS.png", width = 12, height = 8, units = "in", res = 500)
#Mfuzz::mfuzz.plot2(eset = eset_PP.fs, cl = clusters_PP, mfrow = c(2,(c_opt/2)), x11 = FALSE, centre = TRUE, centre.col = "black", time.labels = c("16cell","blastula","gastrula","prototroch","swimming","1week"), 
#                   min.mem = 1.01, centre.lwd = 4, xlab = "Developmental Stage", ylab = "Gene Expression (Change)")
#dev.off()

# MFuzz clusters with all genes
#png(filename = "CLUSTERS_heatmap.png", width = 12, height = 8, units = "in", res = 500)
#Mfuzz::mfuzz.plot2(eset = eset_PP.fs, cl = clusters_PP, mfrow = c(2,(c_opt/2)), x11 = FALSE, centre = TRUE, centre.col = "black", time.labels = c("16cell","blastula","gastrula","prototroch","swimming","1week"), 
#                   min.mem = 0.5, centre.lwd = 4, xlab = "Developmental Stage", ylab = "Gene Expression (Change)")
#dev.off()



