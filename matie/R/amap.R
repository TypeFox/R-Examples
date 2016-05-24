# agram.R
# Author: Ben Murrell, Dan Murrell and Hugh Murrell
# Adapted from Moignard et al, 

# clusters variables in a datafrome using dist = 1-A 
# Plots a heatmap with cluster dendrogram attached.

amap <- function (dataSet, 
                  palette=colorRampPalette(c("blue","white", "red"), space = "rgb"),
                  corAdjusted=FALSE, 
                  method="spearman", 
                  title="Association Clustering",
                  ...) {
  
  # require(gplots)
  # require(cba)
  
  if (corAdjusted) {
    matrix <- as.matrix(dataSet) 
    dC <- cor(matrix, method = method, use="pairwise.complete.obs")
  }
  
  dA <- tap(dataSet) 
  dA <- as.matrix(dA)
  diag(dA) <- 1
  # print(dim(dA))
  
  if (corAdjusted) {
    dP <- dA * sign(dC)
  } else {
    dP <- dA
  }
  
  pd <- as.dist(1-dP) 
  phc <- hclust(pd, method = "average", members=NULL)
  pco <- order.optimal(pd, phc$merge)
  pho <- phc
  pho$merge <- pco$merge
  pho$order <- pco$order
  heatmap.2(dP, Rowv = as.dendrogram(pho), Colv = as.dendrogram(pho), 
            dendrogram="row", trace=c("none"), 
            margins=c(5,5), key = TRUE, 
            symkey=corAdjusted, keysize=2, 
            density.info=c("none"), 
            main=title, 
            symbreaks=TRUE, col=palette)
  
  invisible(NULL)
     
}



# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


