#' spectral clustering
#'
#' @param g sg object. Should be weighted (with weight.sg-function)
#' @param m levels to consider
#' @param K number of assumed clusters
#'
#' @importFrom stats dist kmeans runif
#' @export

spectral.sg<-function(g, m=2, K=3) {

  if(!is(g, "sg")) stop("g not of class sg.")
  if(is.null(g$weights)) stop("No weights in x. Run weight.sg-function.")

  W <- sg2wadj(g)$matrix
  G <- diag(rowSums(W))
  L <- G-W # Laplacian

  E <- eigen(L)
  l <- E$values
  chosen <- order(l)[1:m+1] # drop first 0 value
  v <- E$vectors[,chosen]
  labels<-kmeans(v, K)$cluster
  #'
  note<-paste("K-means from weighted", g$type, "(", g$par,") note:", g$note)
  #
  x <- list(id=labels,
       sgc=as.sgc(split(1:g$N,labels), type="spectral clustering", pars=list(m=m, K=K), note=note),
       v=v, l=l,
       N = g$N)
  class(x) <- "sgspectral"
  x
}

#' plot spectral clustering results
#'
#' @param x spectral.sg result
#' @param data point pattern
#' @param ... ignored
#'
#' @export
plot.sgspectral <- function(x, data, ...) {

  if(missing(data)) stop("Can't plot without pattern data.")

  data <- sg_parse_coordinates(data)

  if(ncol(data)>2) stop("plot only for 2d.")

  par(mfrow=c(1,3))

  plot(x$sgc, data, col=x$id, main="Identified clusters", pch=19)

  plot(1:20, sort(x$l)[1:20], main="20 smallest eigenvalues", xlab="Number", ylab="eigenvalue", col="darkgreen", pch=19)
  abline(h=0, col="gray60", lty=2)

  m <- x$sgc$par$m

  plot(NA, NA, xlab="index", ylab="", xlim=c(0, x$N), ylim=c(0, m), yaxt="n", main="Smallest eigenvectors (not in scale)")
  axis(2, at=1:m-0.5, paste("eigen",1:m+1, sep=""), tick=FALSE)

  for(i in 1:m){
    vi<-x$v[,i]-mean(x$v[,i])
    vi<-vi/(1.1*m*max(abs(vi)))+i-0.5
    points(1:x$N, vi, col=x$id)
  }

  abline(h=c(1:(m-1)))
  abline(h=1:m-.5, col="gray50", lty=2)

}
