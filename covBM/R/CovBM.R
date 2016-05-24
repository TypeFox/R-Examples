#' covBM
#'
#' This is a constructor function for the "covBM" class, representing a Brownian motion
#' component in terms of a continuous variable. The object created is a special type of
#' \code{\link[nlme]{corStruct}}.
#' @param value Numeric argument providing starting value for the scale parameter of
#' Brownian motion process relative to residual error variance for optimisation.
#' @param form A one-sided formula of the form ~t|g, where t represents a continuous
#' variable (usually time) and g represents a grouping factor, i.e. with a separate
#' Brownian motion process modelled at each level.
#' @return An object of class "covBM" and inheriting from "corStruct".
#' @export
#' @examples cov1<-covBM(form=~time|group)
covBM <-  function(value=1,form = ~1) {	
	if(value<=0) {stop("Scale parameter for BM process must be positive")}
	value <- log(value)
	attr(value, "formula") <- form
    attr(value, "fixed") <- FALSE
    class(value) <- c("covBM", "corStruct")
    value
}


makecovBM<-function(variate,kappa){
  val<-.C(covBM_C,as.double(variate),as.double(kappa), as.integer(length(variate)), result=double(length(variate)*length(variate)))$result
  val
  }


#' corMatrix.covBM
#'
#' This method generates a scaled covariance matrix (or list of matrices), for a "covBM" "corStruct" object.
#' @param object An object of class \code{\link{covBM}}, inheriting from \code{\link[nlme]{corStruct}}.
#' @param covariate List of covariate vectors, at which values the correlation matrix,
#' or list of correlation matrices, are to be evaluated, as for \code{\link[nlme]{corMatrix.corStruct}}.
#' @param ... Additional arguments (not used by this method).
#' @export
corMatrix.covBM <- function(object, covariate = getCovariate(object), ...)
{
     	cm <- lapply(covariate,makecovBM,kappa=exp(as.vector(object)))
         return(cm)
}


#' coef.covBM
#'
#' This is a method function that extracts the scale coefficient associated with a
#' Brownian motion correlation structure object.
#' @param object An object of class \code{\link{covBM}}, inheriting from \code{\link[nlme]{corStruct}}.
#' @param unconstrained A logical value. If TRUE the coefficients are returned in unconstrained form (as used
#' in the optimization algorithm). If FALSE the coefficients are returned in "natural" form.
#' @param ... Additional arguments (not used by this method).
#' @export
#' @examples cov1<-covBM(form=~time|group)
#' coef(cov1)
coef.covBM <- function(object, unconstrained=TRUE, ...) {
  if (unconstrained) {
    if (attr(object, "fixed")) {
      return(numeric(0))
    } else {
      return(as.vector(object))
    }
  }
  aux <- exp(as.vector(object))
  names(aux)<-"Kappa"
  aux
}
