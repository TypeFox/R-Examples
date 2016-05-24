##' Boosting ordinal Bradley-Terry-Luce models
##' 
##' Boosting procedure for ordinal BTL models
##' 
##' @usage BTLboost(formula, data, objects=NULL, groupVars=NULL,
##'                 selection=c("DEVIANCE","AIC","BIC"),
##'                 mstop=500, nu=1, maxit=1, verbose=TRUE, ...)
##' 
##' @param formula a formula describing the full model.
##' @param data a data frame containing the design matrix for the model 
##' (See also \code{\link{design}} to generate such an design matrix).
##' @param objects (optional) a character vector specifying the objects that 
##' should always be part of the model.
##' @param groupVars (optional) a character vector specifying the subject-specific 
##' covariates, whose subject-object interactions are considered simultaneously 
##' in each boosting step. By default (\code{groupVars=NULL}), all subject-object
##' interactions are cosidered separately in each boosting step.
##' @param selection a character specifying the criterion that is used in each 
##' boosting step to determine the best fitting covariate(s).
##' @param mstop an integer giving the number of boosting iterations.
##' @param nu a double between 0 and 1 defining the step size or shrinkage parameter.
##' @param maxit an integer representing the maximum number of Fisher-scoring iterations (see also \code{\link[VGAM]{vglm.control}}).
##' @param verbose logical indicating if output should be produced for each boosting iteration.
##' @param ... further arguments passed to \code{ordBTL}.
##'
##' @author Giuseppe Casalicchio
##' 
##' @return A List of \itemize{
##' \item \code{BEST} contains estimated parameters of the last boosting iteration
##' \item \code{AIC} a vector of AIC values for each boosting iteration
##' \item \code{BIC} a vector of BIC values for each boosting iteration
##' \item \code{DEVIANCE} a vector that reflects the deviance of each boosting iteration
##' \item \code{PATH} a dataframe containing the coefficient build-up at the end of each boosting iteration
##' \item \code{UPDATED} a vector of strings containing the selected components in each boosting iteration
##' }
##' 
##' @example inst/examples/BTLboost_ex.R
##' @export 
##' @importFrom gtools smartbind
# 
# BTLboost <- function(formula, data, objects=NULL, groupVars=NULL,
#                      selection=c("DEVIANCE","AIC","BIC"),
#                      mstop=500, nu=1, maxit=1, verbose=TRUE, ...){
#   coefFun <- function(model){
#     co <- coefficients(model)
#     if(any(grepl(":",names(co)))){
#       ind <- grepl("GAMMA",gsub(".*:","", names(co)))
#       names(co)[ind] <- paste(gsub(".*:","", names(co)[ind]), ":",
#                               gsub(":.*","", names(co)[ind]), sep="")
#     }
#     return(co)
#   }
#   if(sum(nu > 1 | nu <= 0)==1) stop("nu should be between 0 and 1")
#   if(!all(objects%in%colnames(data))){
#     stop("Not all 'objects' are column names of 'data'")
#   }
#   selection <- match.arg(selection)
#   
#   # model formulas
#   allVars <- attr(terms(formula),"term.labels")
#   if(is.null(objects)) objects <- allVars[allVars%in%colnames(data)]
#   pool <- allVars[!allVars%in%objects]
#   response <- as.character(formula)[2]
#   formObjects <- paste(response, paste(objects, collapse="+"), sep="~")
#   
#   #subjectVars <- groupVars #unique(unlist(strsplit(pool, ":"))[!unlist(strsplit(pool, ":"))%in%objects])
#   if(!is.null(groupVars)){
#     if(!all(groupVars%in%colnames(data))){
#       stop("Not all 'objects' are columns of 'data'")
#     }
#     singlePool <- pool[-grep(paste(groupVars, collapse="|"), pool)]
#     groupPool <- pool[grep(paste(groupVars, collapse="|"), pool)]
#     groupPoolVars <- gsub("GAMMA.*:", "", groupPool)
#     groupPoolList <- vector("list", length(groupVars))
#     for(i in 1:length(groupVars)) 
#       groupPoolList[[i]] <- paste(groupPool[groupPoolVars%in%groupVars[i]], collapse="+")
#     pool <- c(singlePool, unlist(groupPoolList))
#   }
#   formCandidate <- paste(response, paste(pool, "1", sep="-"), sep="~")
#   
#   # preallocation
#   metric <- nparam <- numeric(length(formCandidate))
#   DEVIANCE <- AIC <- BIC <- selected <- rep(NA, length(mstop))
#   
#   ### Step 1: get initial values (intercept model)
#   formula <- as.formula(paste(response,"~", "1"))
#   intercept <- ordBTL(formula, data=data, 
#                       x.arg = FALSE, y.arg = FALSE, qr.arg = FALSE, ...)
#   f <- intercept@predictors
#   param <- coefFun(intercept)
#   
#   # Start loop for m=1,..., mstop 
#   for(j in 1:mstop){
#     base <- coeflist <- vector("list", length=length(formCandidate))
#     ### Step 2A: object predictors
#     formula <- as.formula(formObjects)
#     obligatory <- ordBTL(formula, data=data, offset=f, etastart=f, maxit=maxit, 
#                          x.arg = FALSE, y.arg = FALSE, qr.arg = FALSE, ...)
#     # update step
#     f <- f + nu*obligatory@predictors
#     coefs <- smartbind(param, nu*coefFun(obligatory))
#     param <- colSums(coefs, na.rm=TRUE)
#     
#     ### Step 2B: candidate predictors
#     for(i in 1:length(formCandidate)) {
#       formula <- as.formula(formCandidate[i])
#       suppressWarnings({base[[i]] <- ordBTL(formula, data=data, offset=f, etastart=f, maxit=maxit, 
#                                             x.arg = FALSE, y.arg = FALSE, qr.arg = FALSE, ...)})
#     }
#     
#     # number of parameters  & metric
#     nparam <- #lapply(as.list(pool), function(X) strsplit(X, "[+]")[[1]])
#       unlist(lapply(base, function(X){length(unique(c(names(coefFun(X)), names(param))))}))
#     
#     coeflist <- lapply(base, function(X) colSums(smartbind(param, nu*coefFun(X)), na.rm=TRUE))
#     
#     dev <- sapply(base, function(X) deviance(X))
#     aic <- dev + 2*nparam
#     bic <- dev + log(nobs(intercept))*nparam
#     
#     if(selection=="AIC"){
#       metric <- aic
#     } else{
#       if(selection=="BIC"){
#         metric <- bic
#       }else{
#         metric <- dev
#       }
#     }
#     
#     # choost best model w.r.t. chosen criteria
#     best <- which.min(metric) 
#     
#     # update step
#     coefs <- smartbind(param, nu*coefFun(base[[best]]))  
#     param <- colSums(coefs, na.rm=TRUE)
#     if(j==1) parampath <- param else parampath <- smartbind(parampath, param)
#     f <- f + nu*base[[best]]@predictors
#     DEVIANCE[j] <- dev[best]
#     BIC[j] <- bic[best]
#     AIC[j] <- aic[best]
#     selected[j] <-  pool[best]
#     
#     if(verbose==TRUE){
#       cat(sprintf("% 5.0f", j),
#           paste(selection, ":", sep=""), sprintf("%.3f", metric[best]),
#           "| vars:", sprintf("% 3.0f", length(param)),
#           "| updated:", selected[j],
#           "\n")
#     }
#   }
#   ret <- list(BEST=parampath, AIC=AIC, BIC=BIC, DEVIANCE=DEVIANCE, 
#               PATH=parampath, UPDATED=selected)
#   if(mstop>1) ret[["BEST"]] <- ret[["BEST"]][which.min(ret[[selection]]),!is.na(ret[["BEST"]][which.min(ret[[selection]]),])]
#   return(ret)
# }

BTLboost <- function(formula, data, objects=NULL, groupVars=NULL,
                     selection=c("DEVIANCE","AIC","BIC"),
                     mstop=500, nu=1, maxit=1, verbose=TRUE, ...){
  coefFun <- function(model){
    co <- coefficients(model)
    if(any(grepl(":",names(co)))){
      ind <- grepl("GAMMA",gsub(".*:","", names(co)))
      names(co)[ind] <- paste(gsub(".*:","", names(co)[ind]), ":",
                              gsub(":.*","", names(co)[ind]), sep="")
    }
    return(co)
  }
  logL <- function(formula, data, coef){
    resp <- as.character(formula)[2]
    if(resp%in%colnames(data)) Y <- data[,resp] else
      stop(paste("'", resp, "'", " not found in columns of 'data'."))
      
    nthres <- length(unique(Y))-1
    X <- model.matrix(formula, data)
    coefNoInt <- coef[!grepl("Intercept", names(coef))]
    
    int <- (coef[grepl("Intercept", names(coef))])
    if(length(int)==0) int <- 0
    if(int[1]>0) int <- -1*int
    
    if(nthres>2*length(int)){ # zero thres missing
      intNew <- c(int,0,-int)
    }else{
      intNew <- c(int, -int)
    }
    int <- intNew
    if(length(coefNoInt)==1){
      xbeta <- X[,names(coefNoInt)]*coefNoInt 
    } else{
      xbeta <- X[,names(coefNoInt)]%*%coefNoInt 
    }
    tmp <- sum((Y==1)*log(exp(xbeta+int[1])/(1+exp(xbeta+int[1]))), na.rm=TRUE)
    for(i in 2:(nthres)){
      tmp <- tmp +sum((Y==i)*(log( exp(xbeta+int[i])/(1+exp(xbeta+int[i])) - exp(xbeta+int[i-1])/(1+exp(xbeta+int[i-1])) )))
    }
    i <- (nthres)+1
    tmp <- tmp + sum((Y==i)*log(1- exp(xbeta+int[i-1])/(1+exp(xbeta+int[i-1]))), na.rm=TRUE)
    tmp
  }
  if(sum(nu > 1 | nu <= 0)==1) stop("nu sh")
  if(!all(objects%in%colnames(data))){
    stop("Not all 'objects' are column names of 'data'")
  }
  args <- list(...)
  if("family"%in%names(args) && args$family!="cumulative" && nu!=1)
    stop("Currently only cumulative logit models with smaller stepsize available")
  if("family.control"%in%names(args) && "link"%in%names(args$family.control) && args$family.control$link!="logit" && nu!=1)
    stop("Currently only cumulative logit models with smaller stepsize available")
  
  selection <- match.arg(selection)
  
  # model formulas
  allVars <- attr(terms(formula),"term.labels")
  if(is.null(objects)) objects <- allVars[allVars%in%colnames(data)]
  pool <- allVars[!allVars%in%objects]
  response <- as.character(formula)[2]
  formObjects <- paste(response, paste(objects, collapse="+"), sep="~")
  
  #subjectVars <- groupVars #unique(unlist(strsplit(pool, ":"))[!unlist(strsplit(pool, ":"))%in%objects])
  if(!is.null(groupVars)){
    if(!all(groupVars%in%colnames(data))){
      stop("Not all 'objects' are columns of 'data'")
    }
    singlePool <- pool[-grep(paste(groupVars, collapse="|"), pool)]
    groupPool <- pool[grep(paste(groupVars, collapse="|"), pool)]
    groupPoolVars <- gsub("GAMMA.*:", "", groupPool)
    groupPoolList <- vector("list", length(groupVars))
    for(i in 1:length(groupVars)) 
      groupPoolList[[i]] <- paste(groupPool[groupPoolVars%in%groupVars[i]], collapse="+")
    pool <- c(singlePool, unlist(groupPoolList))
  }
  formCandidate <- paste(response, paste(pool, "1", sep="-"), sep="~")
  
  # preallocation
  metric <- nparam <- numeric(length(formCandidate))
  DEVIANCE <- AIC <- BIC <- selected <- rep(NA, length(mstop))
  
  ### Step 1: get initial values (intercept model)
  formula <- as.formula(paste(response,"~", "1"))
  intercept <- ordBTL(formula, data=data, 
                      x.arg = FALSE, y.arg = FALSE, qr.arg = FALSE, ...)
  f <- intercept@predictors
  param <- coefFun(intercept)
  
  # Start loop for m=1,..., mstop 
  for(j in 1:mstop){
    base <- coeflist <- vector("list", length=length(formCandidate))
    ### Step 2A: object predictors
    formula <- as.formula(formObjects)
    obligatory <- ordBTL(formula, data=data, offset=f, etastart=f, maxit=maxit, 
                         x.arg = FALSE, y.arg = FALSE, qr.arg = FALSE, ...)
    # update step
    f <- f + nu*obligatory@predictors
    coefs <- smartbind(param, nu*coefFun(obligatory))
    param <- colSums(coefs, na.rm=TRUE)
    
    ### Step 2B: candidate predictors
    for(i in 1:length(formCandidate)) {
      formula <- as.formula(formCandidate[i])
      suppressWarnings({base[[i]] <- ordBTL(formula, data=data, offset=f, etastart=f, maxit=maxit, 
                                            x.arg = FALSE, y.arg = FALSE, qr.arg = FALSE, ...)})
    }
    
    # number of parameters  & metric
    nparam <- #lapply(as.list(pool), function(X) strsplit(X, "[+]")[[1]])
      unlist(lapply(base, function(X){length(unique(c(names(coefFun(X)), names(param))))}))
    
    coeflist <- lapply(base, function(X) colSums(smartbind(param, nu*coefFun(X)), na.rm=TRUE))
    if(nu!=1){
      dev <- sapply(coeflist, function(X) -2*logL(formula, data, X)) #sapply(base, function(X) deviance(X)) #
    }else{
      dev <- sapply(base, function(X) deviance(X))
    }
    aic <- dev + 2*nparam
    bic <- dev + log(nobs(intercept))*nparam
    
    if(selection=="AIC"){
      metric <- aic
    } else{
      if(selection=="BIC"){
        metric <- bic
      }else{
        metric <- dev
      }
    }
    
    # choost best model w.r.t. chosen criteria
    best <- which.min(metric) 
    
    # update step
    coefs <- smartbind(param, nu*coefFun(base[[best]]))  
    param <- colSums(coefs, na.rm=TRUE)
    if(j==1) parampath <- param else parampath <- smartbind(parampath, param)
    f <- f + nu*base[[best]]@predictors
    DEVIANCE[j] <- dev[best]
    BIC[j] <- bic[best]
    AIC[j] <- aic[best]
    selected[j] <-  pool[best]
    
    if(verbose==TRUE){
      cat(sprintf("% 5.0f", j),
          paste(selection, ":", sep=""), sprintf("%.3f", metric[best]),
          "| vars:", sprintf("% 3.0f", length(param)),
          "| updated:", selected[j],
          "\n")
    }
  }
  ret <- list(BEST=parampath, AIC=AIC, BIC=BIC, DEVIANCE=DEVIANCE, 
              PATH=parampath, UPDATED=selected)
  if(mstop>1) ret[["BEST"]] <- ret[["BEST"]][which.min(ret[[selection]]),!is.na(ret[["BEST"]][which.min(ret[[selection]]),])]
  return(ret)
}
