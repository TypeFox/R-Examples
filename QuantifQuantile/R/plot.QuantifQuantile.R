#' @title Plot of estimated conditional quantiles using optimal quantization.
#' @name plot.QuantifQuantile
#' @description This function plots the estimated conditional quantiles by default. 
#' It can also illustrate our data driven selection criterion 
#' for \code{N} by providing the plot of the bootstrap estimated values of
#' integrated squared error ISE(N) versus \code{N}. 
#' @param x An object of class \code{QuantifQuantile}, which is the result of 
#' \code{\link{QuantifQuantile}} or \code{\link{QuantifQuantile.d2}}.
#' @param col.plot Vector of size \code{length(x$alpha)+1}. The first entry 
#' corresponds to the color of the data points while the other colors are for 
#' the conditional quantiles curves, points or surfaces.
#' @param ise Whether it plots the ISE curves in addition to the estimated 
#' quantile curves (if \code{ise=TRUE}, two different plots).
#' @param \dots Arguments to be passed to \code{\link{par}}.
#' 
#' @details If \code{X} is univariate, the graph is two-dimensional and if 
#' \code{X} is bivariate, it provides a 3D-graph using the \code{\link{rgl}} 
#' package. When only one value for \code{x} is considered, estimated 
#' conditional quantiles are plotted as points. When \code{x} is a grid of 
#' values, they are plotted as curves if \code{d}=1 and surfaces if \code{d}=2.
#' 
#' When \code{ise=TRUE}, the first plot allows to adapt the choice of the grid for \code{N}, 
#' called \code{testN}. For example, if the curve is decreasing with \code{N}, it 
#' indicates that the values in \code{testN} are too small and the optimal 
#' \code{N} is larger. 

#' @references Charlier, I. and Paindaveine, D. and Saracco, J.,
#' \emph{Conditional quantile estimation through optimal quantization}, 
#' Journal of Statistical Planning and Inference, 2015 (156), 14-30.
#' @references Charlier, I. and Paindaveine, D. and Saracco, J.,
#' \emph{Conditional quantile estimator based on optimal 
#' quantization: from theory to practice}, Submitted.

#' @seealso \code{\link{QuantifQuantile}}, \code{\link{QuantifQuantile.d2}} and
#'\code{\link{QuantifQuantile.d}}
#' 
#' @examples
#' #for a univariate X
#' set.seed(644972)
#' n <- 300
#' X <- runif(300,-2,2)
#' Y <- X^2+rnorm(n)
#' res <- QuantifQuantile(X,Y,testN=seq(10,25,by=5))
#' plot(res,ise=TRUE)
#' 
#' \dontrun{
#' set.seed(92536)
#' n <- 300
#' X <- runif(300,-2,2)
#' Y <- X^2+rnorm(n)
#' res <- QuantifQuantile(X,Y,testN=seq(10,25,by=5),x=1)
#' plot(res,ise=TRUE)
#' 
#' 
#' #for a bivariate X
#' #(a few seconds to execute)
#' set.seed(253664)
#' d <- 2
#' n <- 1000
#' X<-matrix(runif(d*n,-2,2),nr=d)
#' Y<-apply(X^2,2,sum)+rnorm(n)
#' res <- QuantifQuantile.d2(X,Y,testN=seq(80,130,by=10),B=20,tildeB=15)
#' plot(res,ise=TRUE)
#' 
#' set.seed(193854)
#' d <- 2
#' n <- 1000
#' X<-matrix(runif(d*n,-2,2),nr=d)
#' Y<-apply(X^2,2,sum)+rnorm(n)
#' res <- QuantifQuantile.d2(X,Y,testN=seq(110,140,by=10),x=as.matrix(c(1,0)),
#' B=30,tildeB=20)
#' plot(res,ise=TRUE)
#' }
#' 
#' @importFrom rgl plot3d
#' @importFrom rgl points3d
#' @importFrom rgl surface3d
#' @importFrom grDevices dev.new
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @export plot.QuantifQuantile
#' @S3method plot QuantifQuantile

plot.QuantifQuantile <- function(x, col.plot = c(1:(length(x$alpha) + 1)), ise=FALSE,...) {
    stopifnot(class(x)=="QuantifQuantile")
    if(ise==FALSE){
      if (is.vector(x$X)) {
        plot(x$X, x$Y, col = col.plot[1], cex = 0.7,...)
        if (length(x$x) > 1) {
          for (j in 1:length(x$alpha)) {
            lines(x$x, x$hatq_opt[j, ], col = col.plot[j + 
                                                         1])
          }
        }
        if (length(x$x) == 1) {
          for (j in 1:length(x$alpha)) {
            points(x$x, x$hatq_opt[j, ], col = col.plot[j + 
                                                          1], pch = 20, cex = 1)
          }
        }
      } else {
        d <- nrow(x$X)
        if (d == 1) {
          plot(x$X, x$Y, col = col.plot[1], cex = 0.7, ...)
          if (length(x$x) > 1) {
            for (j in 1:length(x$alpha)) {
              lines(x$x, x$hatq_opt[ j,], col = col.plot[j + 
                                                           1])
            }
          }
          if (length(x$x) == 1) {
            for (j in 1:length(x$alpha)) {
              points(x$x, x$hatq_opt[ j, ], col = col.plot[j + 
                                                             1], pch = 20, cex = 1)
            }
          }
        }
        if(d ==2) {
          plot3d(x$X[1, ], x$X[2, ], x$Y, col = col.plot[1], ...)
          hatq_matrix <- array(0, dim = c(sqrt(dim(x$x)[2]), sqrt(dim(x$x)[2]), 
                                          length(x$alpha)))
          for (i in 1:length(x$alpha)) {
            hatq_matrix[, , i] <- matrix(x$hatq_opt[ i,], ncol = dim(x$x)[2])
          }
          if (length(x$x)/d > 1) {
            for (j in 1:length(x$alpha)) {
              surface3d(unique(x$x[1, ]), unique(x$x[2, ]), 
                        hatq_matrix[, , j], col = col.plot[j + 1])
            }
          }
          if (length(x$x)/d == 1) {
            for (j in 1:length(x$alpha)) {
              points3d(x$x[1, ], x$x[2, ], hatq_matrix[, , 
                                                       j], col = col.plot[j + 1], size = 5)
            }
          }
        }
        if (d > 2) {
          stop("No graphical illustration available when d>2")
        }
      }  
    }else{
      par(mar=c(5, 4.7, 4, 2) + 0.1)
      if(length(x$N_opt)==1){
        hatISEmean_N <- apply(x$hatISE_N, 2, mean)
        plot(x$testN, hatISEmean_N, type = "l", xlab="N",ylab=expression(hat(ISE)(N)))
        abline(v = x$N_opt, col = 2)
      }else{
        plot(x$testN, x$hatISE_N[1,], type = "l", ylim = c(min(x$hatISE_N),max(x$hatISE_N)),col=1,lty=1,
             xlab="N",ylab=expression(hat(ISE)[alpha](N)))
        j=which(x$N_opt[1]==x$testN)
        lines(c(x$N_opt[1],x$N_opt[1]),c(-1,x$hatISE_N[1,j]),
              col=2)
        if(length(x$alpha)>= 1){
            for(i in 2:length(x$alpha)){
              lines(x$testN, x$hatISE_N[i,], type = "l",col=1,lty=i,lwd=1)
              j=which(x$N_opt[i]==x$testN)
              lines(c(x$N_opt[i],x$N_opt[i]),c(-1,x$hatISE_N[i,j]),
                    col=2)
            }
            }
        legend("top",legend = x$alpha,lty=c(1:length(x$alpha)))
          }
      if (is.vector(x$X)) {
        dev.new()
        plot(x$X, x$Y, col = col.plot[1], cex = 0.7,...)
        if (length(x$x) > 1) {
          for (j in 1:length(x$alpha)) {
            lines(x$x, x$hatq_opt[j, ], col = col.plot[j + 
                                                         1])
          }
        }
        if (length(x$x) == 1) {
          for (j in 1:length(x$alpha)) {
            points(x$x, x$hatq_opt[j, ], col = col.plot[j + 
                                                          1], pch = 20, cex = 1)
          }
        }
      } else {
        d <- nrow(x$X)
        if (d == 1) {
          dev.new()
          plot(x$X, x$Y, col = col.plot[1], cex = 0.7, ...)
          if (length(x$x) > 1) {
            for (j in 1:length(x$alpha)) {
              lines(x$x, x$hatq_opt[ j,], col = col.plot[j + 
                                                           1])
            }
          }
          if (length(x$x) == 1) {
            for (j in 1:length(x$alpha)) {
              points(x$x, x$hatq_opt[ j, ], col = col.plot[j + 
                                                             1], pch = 20, cex = 1)
            }
          }
        }
        if(d ==2) {
          plot3d(x$X[1, ], x$X[2, ], x$Y, col = col.plot[1], ...)
          hatq_matrix <- array(0, dim = c(sqrt(dim(x$x)[2]), sqrt(dim(x$x)[2]), 
                                          length(x$alpha)))
          for (i in 1:length(x$alpha)) {
            hatq_matrix[, , i] <- matrix(x$hatq_opt[ i,], ncol = dim(x$x)[2])
          }
          if (length(x$x)/d > 1) {
            for (j in 1:length(x$alpha)) {
              surface3d(unique(x$x[1, ]), unique(x$x[2, ]), 
                        hatq_matrix[, , j], col = col.plot[j + 1])
            }
          }
          if (length(x$x)/d == 1) {
            for (j in 1:length(x$alpha)) {
              points3d(x$x[1, ], x$x[2, ], hatq_matrix[, , 
                                                       j], col = col.plot[j + 1], size = 5)
            }
          }
        }
        if (d > 2) {
          stop("No graphical illustration available when d>2")
        }
      } 
        }   
} 
