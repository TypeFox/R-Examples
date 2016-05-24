#' Plot matrix
#'
#' This function plots a data matrix
#'
#' @param y matrix to be ploted
#' @param standardize.rows standardize matrix rows for plot
#' @param reorder.rows reorder matrix rows based on pattern
#' @param reorder.cols reorder matrix cols based on pattern
#' @param high.contrast apply transformation to matrix to increase contrast
#' @importFrom lattice levelplot
#' @importFrom grDevices colorRampPalette
#' @importFrom stats median
#' @importFrom stats quantile
#' @export

plot_matrix <- function(y, standardize.rows = TRUE, reorder.rows = TRUE, reorder.cols = TRUE, high.contrast = TRUE) {

  y <- as.matrix(y)
  xl <- "Microarrays" 
  yl <- "Probes"
  
  if(standardize.rows) {
    y <- t(apply(y,1,scale)) # standardize the rows of y
  } 
  
  if(reorder.rows) {
    medrow <- apply(y,1,median) # median of the rows of y
    y <- y[order(medrow),] # order the rows of y 
  } 
  
  if(reorder.cols) {
    medcol <- apply(y,2,median) # median of the columns of y
    y <- y[,order(medcol)]  # order the columns of y
  }
  
  if(high.contrast) {
    y <- (abs(y)^(1/3))*sign(y)
  } # higher contrast
  
  mi <- as.numeric(quantile(y,0.001))
  ma <- as.numeric(quantile(y,0.999))
  nr <- nrow(y)
  nc <- ncol(y)

  if(nr==nc) { # Image of Correlation and Covariance matrices
   xl <- ""
   yl <- ""
   mi <- ifelse(ma==1, -1, ma)
  }

  spr <- ifelse(nr<=5, 1, round(nr/10))
  spc <- ifelse(nc<=5, 1, round(nc/10))
  sc <- list(
    x = list(at=c(seq(1,nc, spc),nc), labels=c(seq(1,nc,spc),nc)), 
    y = list(draw=FALSE)
  )
  col.l <- colorRampPalette(c('blue','white','red'))
  cbar <- seq(mi,ma,length.out=100)
  
  levelplot(t(y[nr:1,]),col.regions=col.l,xlab=xl,ylab=yl,scales=sc,at=cbar,aspect="fill")

}