#' covFracBM
#'
#' This is a constructor function for the "covFracBM" class, representing a fractional Brownian
#' motion component in terms of a continuous variable. The object created is a special type of
#' \code{\link[nlme]{corStruct}}. 
#' @param value Vector of length 2 providing starting values for optimisation of the scale
#' parameter of fractional Brownian motion process relative to residual error variance and
#' the Hurst parameter, respectively.
#' @param form A one-sided formula of the form ~t|g, where t represents a continuous
#' variable (usually time) and g represents a grouping factor, i.e. with a separate
#' fractional Brownian motion process modelled at each level.
#' @return An object of class "covFracBM" and inheriting from "corStruct".
#' @export
#' @examples cov2<-covFracBM(form=~time|group)
covFracBM <-  function(value=c(1,0.5),form = ~1) {
	if(value[1]<=0) {stop("Scale parameter for fBM process must be positive")}
	if(value[2]<=0 | value[2]>=1) {stop("Hurst parameter for fBM process must be between 0 and 1")}	
	value[1] <- log(value[1])
	value[2]<-logit(value[2])
    attr(value, "formula") <- form
    attr(value, "fixed") <- FALSE
    class(value) <- c("covFracBM", "corStruct")
    value
}


logit <- function(x) { log((x)/(1-x)) }


logit_inv <- function(x) {
  ex <- exp(x)
  ex/(1+ex)
}


makecovFracBM<-function(variate, kappa, Hurst){
  val<-.C("covFracBM_C",as.double(variate),as.double(kappa), as.double(Hurst), as.integer(length(variate)), result=double(length(variate)*length(variate)))[["result"]]
  val
  }

#' corMatrix.covFracBM
#'
#' This method generates a scaled covariance matrix (or list of matrices), for a "covFracBM" "corStruct" object.
#' @param object An object of class \code{\link{covFracBM}}, inheriting from \code{\link[nlme]{corStruct}}.
#' @param covariate List of covariate vectors, at which values the correlation matrix,
#' or list of correlation matrices, are to be evaluated, as for \code{\link[nlme]{corMatrix.corStruct}}.
#' @param ... Additional arguments (not used by this method).
#' @export
corMatrix.covFracBM <- function(object, covariate = getCovariate(object), ...)
{
     	cm <- lapply(covariate,makecovFracBM,kappa=exp(as.vector(object[1])),Hurst=logit_inv(as.vector(object[2])))
         return(cm)
}


#' coef.covFracBM
#'
#' This is a method function that extracts the scale coefficient and Hurst parameter
#' associated with a fractional Brownian motion correlation structure object.
#' @param object An object of class \code{\link{covFracBM}}, inheriting from \code{\link[nlme]{corStruct}}.
#' @param unconstrained A logical value. If TRUE the coefficients are returned in unconstrained form (as used
#' in the optimization algorithm). If FALSE the coefficients are returned in "natural" form.
#' @param ... Additional arguments (not used by this method).
#' @export
#' @examples cov2<-covFracBM(form=~time|group)
#' coef(cov2)
coef.covFracBM <- function(object, unconstrained=TRUE, ...) {
  if (unconstrained) {
    if (attr(object, "fixed")) {
      return(numeric(0))
    } else {
      return(as.vector(object))
    }
  }
  aux<-numeric(2)
  aux[1] <- exp(as.vector(object[1]))
  aux[2]<-logit_inv(as.vector(object[2]))
  names(aux)<-c("Kappa","Hurst index")
  aux
}

