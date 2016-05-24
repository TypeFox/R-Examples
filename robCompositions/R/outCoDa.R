#' Outlier detection for compositional data
#' 
#' Outlier detection for compositional data using standard and robust
#' statistical methods.
#' 
#' The outlier detection procedure is based on (robust) Mahalanobis distances
#' after a isometric logratio transformation of the data.  Observations with
#' squared Mahalanobis distance greater equal a certain quantile of the
#' Chi-squared distribution are marked as outliers.
#' 
#' If method \dQuote{robust} is chosen, the outlier detection is based on the
#' homogeneous majority of the compositional data set.  If method
#' \dQuote{standard} is used, standard measures of location and scatter are
#' applied during the outlier detection procedure.
#' 
#' plot method: the Mahalanobis distance are plotted against the index.
#' The dashed line indicates the (1 - alpha) quantile of the Chi-squared
#' distribution. Observations with Mahalanobis distance greater than this
#' quantile could be considered as compositional outliers.
#' 
#' @aliases outCoDa print.outCoDa plot.outCoDa
#' @param x compositional data
#' @param quantile quantile, corresponding to a significance level, is used as
#' a cut-off value for outlier identification: observations with larger
#' (squared) robust Mahalanobis distance are considered as potential outliers.
#' @param method either \dQuote{robust} (default) or \dQuote{standard}
#' @param h the size of the subsets for the robust covariance estimation
#' according the MCD-estimator for which the determinant is minimized (the
#' default is (n+p+1)/2).
#' @param coda if TRUE, data transformed to coordinate representation before outlier detection. 
#' @param y unused second plot argument for the plot method
#' @param ... additional parameters for print and plot method passed through
#' @return \item{mahalDist }{resulting Mahalanobis distance} \item{limit
#' }{quantile of the Chi-squared distribution} \item{outlierIndex }{logical
#' vector indicating outliers and non-outliers} \item{method }{method used}
#' @note It is highly recommended to use the robust version of the procedure.
#' @author Matthias Templ, Karel Hron
#' @seealso \code{\link{isomLR}}
#' @references Egozcue J.J., V. Pawlowsky-Glahn, G. Mateu-Figueras and C.
#' Barcel'o-Vidal (2003) Isometric logratio transformations for compositional
#' data analysis. \emph{Mathematical Geology}, \bold{35}(3) 279-300. \
#' 
#' Filzmoser, P., and Hron, K. (2008) Outlier detection for compositional data
#' using robust methods. \emph{Math. Geosciences}, \bold{40} 233-248.\
#' 
#' Rousseeuw, P.J., Van Driessen, K. (1999) A fast algorithm for the minimum
#' covariance determinant estimator.  \emph{Technometrics}, \bold{41} 212-223.
#' @keywords multivariate
#' @export
#' @examples
#' 
#' data(expenditures)
#' oD <- outCoDa(expenditures)
#' oD
#' ## providing a function:
#' oD <- outCoDa(expenditures, coda = log)
#' 
outCoDa <- function(x, quantile=0.975, method="robust", 
                    h=1/2, coda = TRUE){
	if(dim(x)[2] < 2) stop("need data with at least 2 variables")
	
	covEst <- function(x, type) {
		standard <- function(x){
				list(mean=colMeans(x, na.rm=TRUE), varmat=cov(x))  
		}
		robust <- function(x){
				v <- robustbase::covMcd(x)
				list(mean=v$center, varmat=v$cov)
		}
		switch(type,
				standard = standard(x),
				robust = robust(x))
	}
	if(!is.logical(coda) & !is.function(coda)){
	  stop("coda must be logical or function")
	}
	if(!is.logical(coda)){
	  x <- coda(x)
	}	else if (coda){ 
	  x <- isomLR(x)
	}  
	cv <- covEst(x, "robust")
	cvc <- covEst(x, "standard")
	dM <- sqrt(mahalanobis(x, center=cv$mean, cov=cv$varmat))
	dMc <- sqrt(mahalanobis(x, center=cvc$mean, cov=cvc$varmat))
	limit <- sqrt(qchisq(p=quantile, df=ncol(x)-1))
	res <- list(mahalDist = dM, limit = limit, 
			    outlierIndex = dM > limit, method=method, 
			    om2 = dMc > limit,
			    m2=dMc, coda = coda)
	class(res) <- "outCoDa"
    invisible(res)
}

#' @rdname outCoDa
#' @export
#' @method print outCoDa
print.outCoDa <- function(x, ...){
  cat("\n --------------------\n")	
  print(paste(length(which(x$outlierIndex == TRUE)), "out of", length(x$outlierIndex), "observations are detected as outliers."))
  cat("\n --------------------\n\n")		
}
#' @rdname outCoDa 
#' @export
#' @method plot outCoDa
#' @param which 1 ... MD against index
#' 2 ... distance-distance plot
plot.outCoDa <- function(x, y, ..., which = 1){
  #	plot(1:length(x$mahalDist), x$mahalDist, ylab="Mahalanobis distance", xlab="Index", type="n", ylim=c(0, max(x$mahalDist)) )
  #	points(1:length(x$mahalDist[x$outlierIndex]), x$mahalDist[x$outlierIndex], pch=3, col="red")
  #	points(1:length(x$mahalDist[!x$outlierIndex]), x$mahalDist[!x$outlierIndex])	
  #	abline(h=x$limit, lty=2)
  #	#print(x$mahalDist)
  if(which == 1){
  if(x$method =="standard") yl <- "Mahalanobis distance" else yl <- "Robust Mahalanobis distance" 
    plot(1:length(x$mahalDist), x$mahalDist, ylab = yl,
         xlab = "Index", ylim = c(0, max(x$mahalDist)),
         pch = as.numeric(x$outlierIndex)*2+1)
    abline(h=x$limit, lty=2)
  }
  if(which == 2){
    if(x$method == "standard"){
      return(cat("apply method robust in outCoDa to make a distance-distance plot available"))
    }
    n <- length(x$mahalDist)
    outlier <- ifelse(x$outlierIndex & x$om2, "both", 
                      ifelse(x$outlierIndex & !x$om2, "robust only",
                             ifelse(x$om2 & !x$outlierIndex, "standard only", "no outlier")))
    MD <- NULL
    RMD <- NULL
    limit <- NULL
    name <- NULL
    w <- NULL
    df <- data.frame("RMD" = x$mahalDist,
                     "MD" = x$m2,
                     "limit" = rep(x$limit, n),
                     "outlier" = outlier,
                     "name" = 1:n)
    gg <- ggplot(df, aes(x = MD, y = RMD, colour = outlier )) 
    gg <- gg + geom_point(data=subset(df, RMD <= limit & MD <= limit)) 
    gg <- gg + geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour="grey") 
    gg <- gg + geom_text(data=subset(df, RMD > limit | MD > limit),  aes(MD, RMD, label=name))
    gg <- gg + geom_vline(xintercept = x$limit, size=0.25) + geom_hline(yintercept = x$limit, size=0.25)
    gg <- gg + theme_bw()
    print(gg)
  }
}
