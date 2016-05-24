#' @title Distribution Selection Criteria.
#' 
#' @description A function to calculate the distribution selection criteria 
#' for a list of candidate fits.
#'
#' @rdname DistSelCriteria
#' @name DistSelCriteria
#' @aliases DistSelCriteria
#' 
#' @details When comparing models fitted by maximum likelihood to the same data, the smaller the AIC, BIC or MDL, the better the fit.
#' When comparing models using the log-likelihood criterion, the larger the log-likelihood the better the fit.
#'
#' @note The MDL criterion only works for parameter estimation by numerical maximum likelihood.
#' 
#' @param X Sample obersevations.
#' @param w An optional vector of sample weights.
#' @param candDist A vector of names of candidate distributions.
#' @param criteria A vector of criteria to be calculated.
#'
#' @author Haizhen Wu and A. Jonathan R. Godfrey.
#' 
#' @return An object of class matrix, containing the listed distribution selection criteria for the named distributions.
#' @export DistSelCriteria
#'
#' @examples 
#' Ozone <- airquality$Ozone
#' Ozone <- Ozone[!is.na(Ozone)] # Removing the NA's from Ozone data
#' DistSelCriteria(Ozone, candDist = c("Gamma", "Weibull", "Normal", "Exp"),
#' criteria = c("logLik","AIC","AICc", "BIC"))

#' @name DistSelCriteria
#' @export DistSelCriteria
DistSelCriteria <- function(X,w = rep(1,length(X))/length(X),candDist = c("Beta_ab","Laplace","Normal"), criteria = c("logLik","AIC","AICc","BIC","MDL")){
    
    if(!(all(criteria %in% c("logLik","AIC","AICc","BIC","MDL")))) {return("some criteria unknown")}
    
    w <- w/sum(w)*length(X)   
    est.pars <- lapply(paste0("e",candDist), do.call, args=list(X,w)) 

    
    function.list <- lapply(criteria, getS3method, class = "eDist")
    fn <- function (est.par) (lapply(X=function.list, FUN=do.call, args=list(est.par)))
    CriteriaValues <- sapply(est.pars, fn )
    colnames(CriteriaValues) <- candDist
    rownames(CriteriaValues) <- criteria
    
    return(CriteriaValues)
  }

