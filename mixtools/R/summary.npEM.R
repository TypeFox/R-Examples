summary.npEM <- function(object, ...) {
  normpost <- sweep(object$post, 2, sums <- colSums(object$post), "/")
  n <- NROW(object$data)
  r <- NCOL(object$data)  
  m <- length(object$lambdahat)
  B <- 1
  if (r>1) B <- max(object$blockid)
  lambda <- sums/n
  means <- variances <- NULL
  for(i in 1:B) {    
    if (r>1) {
    	coords <- object$blockid == i
    	xx <- as.vector(object$data[,coords])
    	}
    	else {
    		coords <- 1
    		xx <- as.vector(object$data)
    		}
    sc <- sum(coords)
    M <- V <- NULL
    for (j in 1:m) {
      wts <- rep(normpost[,j]/sc, sc)
      M <- c(M, tmp <- sum(xx*wts))
      V <- c(V, sum((xx-tmp)^2 *wts))
#      cat(M," ")
    }
    means <- rbind(means, M)
    variances <- rbind(variances, V)
  }
  rownames(means) <- rownames(variances) <- paste("block", 1:B)
  colnames(means) <- colnames(variances) <- paste("component", 1:m)
  ans <- list(n=n, m=m, r=r, B=B, blockid=object$blockid,
              means = means, variances=variances)
  class(ans) <- "summary.npEM"
  ans
}

print.summary.npEM <- function(x, digits=3, ...) {	
	if (x$r>1) cat (paste(x$n,"observations,",x$r,"coordinates,",
  		x$m,"components, and",x$B,"blocks.\n\n"))
  	else cat(paste(x$n,"univariate observations, and",
  		x$m,"components.\n\n"))
   		
  cat ("Means (and std. deviations) for each component:\n")
  for(i in 1:x$B) {
  	coords <- 1
  	if (x$r>1) {
  		coords <- x$blockid == i
  		cat(paste("  Block #",i,":  Coordinate", sep=""))
  		cat(ifelse(sum(coords)>1, "s ", " "))
  		   cat(which(coords))
    	cat("\n    ")
  		}   
    cat(paste(signif(x$means[i,],digits), 
              " (", signif(sqrt(x$variances[i,]),digits), ")  ", sep=""))
    cat("\n")
  }
}

