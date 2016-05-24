#' Plot CEC
#' @export
#' @method plot Rcpp_CecModel
#'
#' @title plot
#' 
#' @rdname plot.cec
#' 
#' @description Plot clustering found on 2D plot coloring by cluster.
#' 
#' 
#' @param x CEC model object.
#' @param slice List of dimentions chosen for display since plot is 2D.
#' @param ellipses Outline clusters.
#' @param centers Marks center of every cluster.
#' @param ... other arguments not used by this method.
#' 
#' @param pca Apply PCA or not
#' 
#' @examples
#' \dontrun{
#' plot(cec)
#' plot(cec, slice=c(1,3), ellipses=TRUE)
#' plot(cec, slice=c(1,2,3))
#' plot(cec, ellipses=TRUE, centers=FALSE)
#' plot(cec, pca=TRUE, ellipses=TRUE, centers=FALSE)
#' }
plot.Rcpp_CecModel <- function(x, slice = c(), pca=FALSE, ellipses = FALSE, centers = FALSE, ...) {
  
  d <- x$.getDataset()
  if(pca){
    if(ncol(d) <= 2){
      stop("CEC dataset should have dimension > 2 to use PCA")
    }
    mx <- colMeans(d)
    pca_data <- prcomp(d, scale=FALSE)
    v <- pca_data$rotation
    v <- v[, 1:2]
    d <- pca_data$x
  }
  if (length(slice) == 0) {
    if(pca){
      slice <- c(1,2)
    } else {
      slice <- c(1:(dim(d)[2]))        
    }
    plot(d[,slice], col = (x$clustering + 1), pch=20)
  }
  else if (length(slice) == 1 || length(slice) == 2) {
    plot(d[,slice], col = (x$clustering + 1), pch=20)
  }
  else{
    pairs(d[,slice], col = (x$clustering + 1))
  }
  
  if (ellipses || centers) {
    cen <- x$centers
    n <- length(cen)
    if(pca){
      for (i in 1:n) {
        # t(t(cen[[i]])) creates vector, t(v) is rotation matrix to lower dim subspace
        # t(everything) makes it again a row 
        cen[[i]] <- t(t(v) %*% t(t(cen[[i]])))
      }
    }
    if (ellipses && length(slice) <= 2){
      #library("car")
      cov <- x$covMatrix     
      for (i in 1:n) {
        data <- unlist(cov[i])
        covMat <- matrix(data,ncol=sqrt(length(data)))
        if(pca){
          covMat <- t(v) %*% covMat %*% v
        } else {
          covMat <- covMat[slice,slice]
        }          
        m <-unlist(cen[i][slice])
        eigenValuesAndVectors <- eigen(covMat)
        veE <- eigenValuesAndVectors$vectors
        l <- eigenValuesAndVectors$values
        r <- seq(-pi, pi, by = 0.001)
        len <- length(r)
        Xa <- 2*sqrt(l[1])*cos(r)
        Ya <- 2*sqrt(l[2])*sin(r)
        mm <- c(rep(m[1], len),rep(m[2],len))
        meansMultiply <- matrix(mm, ncol = 2)
        line1 <- cbind(Xa,Ya)
        lineAll <- rbind(line1)
        ddd <- (lineAll%*%t(veE)) + meansMultiply
        points(ddd,col = "black", type = "l", lwd = 2)
        #dataEllipse(d[x$clustering() == (i-1),], plot.points=FALSE, add = TRUE, levels = c(0.9))
      }        
    }
    
    if(centers) {
      mcenters <- do.call(rbind,cen)
      points(mcenters[,slice], col="blue", bg=par("bg"))
    }
  }
}
