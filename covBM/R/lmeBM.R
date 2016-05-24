#' lmeBM
#'
#' This function is a wrapper for \code{\link[nlme]{lme.formula}} that allows
#' Brownian motion, fractional Brownian motion or integrated Ornstein-Uhlenbeck
#' components to be included in linear mixed models, with related parameter
#' estimates and confidence intervals returned in their natural parameterisation.
#' @param fixed This is as specified for \code{\link[nlme]{lme.formula}}.
#' @param random This is as specified for \code{\link[nlme]{lme.formula}}.
#' @param data This is as specified for \code{\link[nlme]{lme.formula}}.
#' @param covariance An optional \code{\link[nlme]{corStruct}} object describing the within-group
#' covariance structure. In addition to those available in \code{nlme},
#' \code{\link{covBM}} can be used to incorporate a Brownian motion component, \code{\link{covFracBM}}
#' can be used to incorporate a fractional Brownian motion component and \code{\link{covIOU}}
#' can be used to incorporate an integrated Ornstein-Uhlenbeck process in relation to 
#' a continuous variable.
#' @param method This is as specified for \code{\link[nlme]{lme.formula}}.
#' @param control This is as specified for \code{\link[nlme]{lme.formula}}.
#' @param keep.data This is as specified for \code{\link[nlme]{lme.formula}}.
#' @return An object of class "lme" representing the linear mixed effects model fit.
#' @export
#' @examples BMmodel<-lmeBM(sqrtcd4~t, data=cd4, random=~t|newpid, covariance=covBM(form=~t|newpid),
#'							method="ML", control=list(opt="link{nlm}"))
lmeBM <-
  function(fixed,
	   data,					
	   random,					
	   covariance = NULL,
	   method = c("REML", "ML"),
	   control = list(),
       keep.data = TRUE)
	{
	origCall<-match.call()
	CallList<-as.list(match.call())[-1]
		
	if(missing(covariance)) stop("lmeBM used without specification of covariance structure")
	if(!missing(control) && !is.list(control)) stop("'control' argument is not provided as a list")
	
	index<-which(names(CallList)=="covariance")
	covObject<-eval(CallList[[index]])
	if(!inherits(covObject, c("covBM", "covFracBM", "covIOU"))) stop("lmeBM used without supported covariance structure")
	names(CallList)[index]<-"correlation"		### covariance argument is fed into 'lme' function as correlation argument, hence renaming here
	
	lme_fit<-do.call(lme, CallList)
	if(!is.null(lme_fit$modelStruct$corStruct) && inherits(lme_fit$modelStruct$corStruct, c("covBM", "covFracBM", "covIOU"))){
		attr(lme_fit$modelStruct$corStruct, "sigma") <-	lme_fit$sigma
		}
		
	if(!is.character(lme_fit$apVar)){lme_fit$apVar<-refactor_apVar(lme_fit$apVar, class(lme_fit$modelStruct$corStruct)[1])}

	lme_fit$call<-origCall

	return(lme_fit)
	}

