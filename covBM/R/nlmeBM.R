#' nlmeBM
#'
#' This function is a wrapper for \code{\link[nlme]{nlme.formula}} that allows
#' Brownian motion, fractional Brownian motion or integrated Ornstein-Uhlenbeck
#' components to be included in non-linear mixed models, with related parameter
#' estimates and confidence intervals returned in their natural parameterisation.
#' @param model This is as specified for \code{\link[nlme]{nlme.formula}}.
#' @param random This is as specified for \code{\link[nlme]{nlme.formula}}.
#' @param start This is as specified for \code{\link[nlme]{nlme.formula}}.
#' @param verbose This is as specified for \code{\link[nlme]{nlme.formula}}.
#' @param control This is as specified for \code{\link[nlme]{nlme.formula}}.
#' @inheritParams lmeBM
#' @return An object of class "nlme" and inheriting from class "lme" representing the non-linear mixed effects model fit.
#' @export
#' @examples data(Milk, package="nlme")
#' Model_fit<- nlmeBM(protein ~ SSasymp(Time, Asym, R0, lrc), data=Milk, 
#'				fixed = Asym + R0 + lrc ~ 1, random = Asym ~ 1|Cow,
#'				covariance=covFracBM(form=~Time|Cow),
#'         		start = c(Asym = 3.5, R0 = 4, lrc = -1))
nlmeBM <-
function(model,
	   data,
	   fixed,
	   random,
	   start,
       covariance = NULL,
 	   method = c("ML", "REML"),
	   control = list(),
	   verbose= FALSE)
	{
	origCall<-match.call()
	CallList<-as.list(match.call())[-1]
		
	if(missing(covariance)) stop("lmeBM used without specification of covariance structure")
	if(!missing(control) && !is.list(control)) stop("'control' argument is not provided as a list")
	
	index<-which(names(CallList)=="covariance")
	covObject<-eval(CallList[[index]])
	if(!inherits(covObject, c("covBM", "covFracBM", "covIOU"))) stop("lmeBM used without supported covariance structure")
	names(CallList)[index]<-"correlation"		### covariance argument is fed into 'lme' function as correlation argument, hence renaming here
	
	nlme_fit<-do.call(nlme, CallList)
	if(!is.null(nlme_fit$modelStruct$corStruct) && inherits(nlme_fit$modelStruct$corStruct, c("covBM", "covFracBM", "covIOU"))){
		attr(nlme_fit$modelStruct$corStruct, "sigma") <-	nlme_fit$sigma
		}
		
	if(!is.character(nlme_fit$apVar)){nlme_fit$apVar<-refactor_apVar(nlme_fit$apVar, class(nlme_fit$modelStruct$corStruct)[1])}

	nlme_fit$call<-origCall

	return(nlme_fit)
	}
