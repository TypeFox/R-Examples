NULL
#'
#' This function makes an inverse Gaussianization procedure besad on PCA iteration ( see \code{\link{inv_GPCA_iteration}}
#'
#' @param x gaussian random variable to transform 
#' @param GPCA_param \code{\link{GPCA-class}} S3 object returned by the function \code{\link{GPCA}}
#' @param extremes see \code{\link{normalizeGaussian_severalstations}}
#' @param type see \code{\link{normalizeGaussian_severalstations}}
#' 
#' @author Emanuele Cordano 
#' 
#' @export 
#' @return the non-Gaussian random variable
#' 
#' @seealso \code{\link{GPCA}},\code{\link{GPCA_iteration}},\code{\link{inv_GPCA_iteration}},\code{\link{inv_GPCA}}
#' @note This function re-iterates the inverse of equation (1) of "PCA Gaussianization for One-Class Remote Sensing Image" by V. Laparra et al.,  \url{http://dx.doi.org/doi/10.1117/12.834011} 
#' @examples 
#' library(RMAWGEN)
#' set.seed(1222)
#' nIterations <- 30
#' N <- 20
#' x <- rexp(N)
#' y <- x+rnorm(N)
#' df <- data.frame(x=x,y=y)
#' 
#' GPCA <- GPCA(df,n=nIterations,extremes=TRUE)
#' 
#' x <- rnorm(N)
#' y <- x+rnorm(N)
#' dfn <- data.frame(x=x,y=y)
#' 
#' GPCAn <- GPCA(dfn,n=nIterations,extremes=TRUE)
#' 
#' df_out <- inv_GPCA(GPCA_param=GPCA,extremes=TRUE)
#' dfn_out <- inv_GPCA(GPCA_param=GPCAn,extremes=TRUE)
#' 
#' 


inv_GPCA <- function (x=NULL,GPCA_param,type=3,extremes=TRUE) {
	
	
	
	n <- length(GPCA_param)-1
	
	if (n<=0) return(x)
	if (is.null(x)) x <- GPCA_param$final_results
	#x <- y ## GUARDARE normalizeGaussian_severalstations 
	y <- normalizeGaussian_severalstations(x,data=GPCA_param[[n]]$x_next,extremes=extremes,inverse=TRUE)  ## ADDED EC 210121206
	for (i in 1:n) {
	
		out <- inv_GPCA_iteration(y,GPCA_param[[n-i+1]],type=type,extremes=extremes)
		y <- out
	}
	

	return(out)

}

