#' Visualize fuzzy clustering results
#' 
#' Takes in a time-course matrix and its clustering results as a cmeans clustering object. Produce a plot to visualize the clustering results.
#' @param Tc a numeric matrix to be clustered. The columns correspond to the time-course and the rows correspond to phosphorylation sites.
#' @param clustObj the clustering of Tc generated from cmeans or kmeans clustering.
#' @param mfrow control the subplots in graphic window.
#' @param cols color palette to be used for plotting. If the color argument remains empty, the default palette is used.
#' @param min.mem phosphorylation sites with membership values below min.mem will not be displayed.
#' @param new.window should a new window be opened for graphics.
#' @param llwd line width. Default is 3.
#' @export
#' 
#' @examples
#' # load the human ES phosphoprotoemics data (Rigbolt et al. Sci Signal. 4(164):rs3, 2011)
#' data(hES)
#' # apply cmeans clustering to partition the data into 11 clusters
#' library(e1071)
#' clustObj <- cmeans(hES, centers=11, iter.max=50, m=1.25)
#' # visualize clustering reuslts
#' fuzzPlot(hES, clustObj, mfrow = c(3,4))
#' 
fuzzPlot <- function (Tc, clustObj, mfrow = c(1, 1), cols, min.mem = 0, new.window = FALSE, llwd=3) {
  clusterindex <- clustObj$cluster
  memship <- clustObj$membership
  memship[memship < min.mem] <- -1
  colorindex <- integer(dim(Tc)[[1]])
  if (missing(cols)) {
    cols <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", 
              "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", 
              "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", 
              "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", 
              "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", 
              "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", 
              "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", 
              "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", 
              "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", 
              "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", 
              "#FF0030", "#FF0018")
  }
  colorseq <- seq(0, 1, length = length(cols))
  
  for (j in 1:max(clusterindex)) {
    tmp <- Tc[clusterindex == j, ]
    tmpmem <- memship[clusterindex == j, j]
    if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0) {
      if (new.window) 
        dev.new()
      par(mfrow = mfrow)
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      plot(x = NA, xlim = c(1 -0.5, dim(Tc)[[2]] +0.5), 
           ylim = c(ymin, ymax), xlab = "Time Course", ylab = "Standardized Profile", 
           main = paste("Cluster", j, "; size=", nrow(tmp)))
    }
    else {
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      plot(x = NA, xlim = c(1, dim(Tc)[[2]]), 
           ylim = c(ymin, ymax), xlab = "Time Course", ylab = "Standardized Profile", 
           main = paste("Cluster", j, "; size=", nrow(tmp)))
    }
    if (!(sum(clusterindex == j) == 0)) {
      for (jj in 1:(length(colorseq) - 1)) {
        tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= 
                     colorseq[jj + 1])
        if (sum(tmpcol) > 0) {
          tmpind <- which(tmpcol)
          for (k in 1:length(tmpind)) {
            lines(tmp[tmpind[k], ], col = cols[jj], lwd=llwd)
          }
        }
      }
    }
  }
}

