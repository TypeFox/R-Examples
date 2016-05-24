#' covIOU
#'
#' This is a constructor function for the "covIOU" class, representing an integrated Ornstein-Uhlenbeck (IOU)
#' process component in terms of a continuous variable. The object created is a special type of
#' \code{\link[nlme]{corStruct}}.
#' @param value Vector of length 2 providing starting values for optimisation of the perturbation
#' parameter of integrated Ornstein-Uhlenbeck process relative to residual error variance and
#' the Alpha parameter, respectively.
#' @param form A one-sided formula of the form ~t|g, where t represents a continuous
#' variable (usually time) and g represents a grouping factor, i.e. with a separate
#' integrated Ornstein-Uhlenbeck process modelled at each level.
#' @return An object of class "covIOU" and inheriting from "corStruct".
#' @export
#' @examples cov3<-covIOU(form=~time|group)
covIOU <-  function(value=c(1,1),form = ~1) {
	if(value[1]<=0) {stop("Perturbation parameter for IOU process must be positive")}
	if(value[2]<=0) {stop("Alpha parameter for IOU process must be positive")}
	value[1] <- log(value[1])
	value[2] <- log(value[2])
    attr(value, "formula") <- form
    attr(value, "fixed") <- FALSE
    class(value) <- c("covIOU", "corStruct")
    value
}

makecovIOU<-function(variate,kappa,alpha){				#covariance matrix for IOU process+error term
	n<-length(variate)
	m<-numeric(n^2)
	rows<-rep(variate,n)
	cols<-rep(variate,rep(n,n,))
	m<- 2*alpha*pmin.int(rows,cols)+exp(-alpha*rows)+exp(-alpha*cols)-1-exp(-alpha*abs(rows-cols))
	(kappa/(2*alpha^3))*matrix(m,n)+diag(n)
}

#' corMatrix.covIOU
#'
#' This method generates a scaled covariance matrix (or list of matrices), for a "covIOU" "corStruct" object.
#' @param object An object of class \code{\link{covIOU}}, inheriting from \code{\link[nlme]{corStruct}}.
#' @param covariate List of covariate vectors, at which values the correlation matrix,
#' or list of correlation matrices, are to be evaluated, as for \code{\link[nlme]{corMatrix.corStruct}}.
#' @param ... Additional arguments (not used by this method).
#' @export
corMatrix.covIOU <- function(object, covariate = getCovariate(object), ...)
{
	cm <- lapply(covariate,makecovIOU,kappa=exp(as.vector(object[1])),alpha=exp(as.vector(object[2])))
      return(cm)
}


#' coef.covIOU
#'
#' This is a method function that extracts the perturbation and Alpha parameters
#' associated with an integrated Ornstein-Uhlenbeck (IOU) process correlation structure object.
#' @param object An object of class \code{\link{covIOU}}, inheriting from \code{\link[nlme]{corStruct}}.
#' @param unconstrained A logical value. If TRUE the coefficients are returned in unconstrained form (as used
#' in the optimization algorithm). If FALSE the coefficients are returned in "natural" form.
#' @param ... Additional arguments (not used by this method).
#' @export
#' @examples cov3<-covIOU(form=~time|group)
#' coef(cov3)
coef.covIOU <- function(object, unconstrained=TRUE, ...) {
  if (unconstrained) {
    if (attr(object, "fixed")) {
      return(numeric(0))
    } else {
      return(as.vector(object))
    }
  }
  aux <- exp(as.vector(object))
  names(aux)<-c("Kappa","Alpha")
  aux
}
