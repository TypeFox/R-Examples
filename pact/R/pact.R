#########################################################################
#		                                                                    #
# R functions for the Predictive Analysis of Clinical Trials R package  #
# These functions handle the low dimensional as well as the high        #
# dimensional case for survival or a binary response variable.          #                                  #									#
#                                                                       #
# Authors: Richard Simon and Jyothi Subramanian	                        #
#			                                                                  #           
#      	                                               		              #
# Nov-2012		                                                          #
#   o Initial programs      	                                          #
# Feb-21-2013                                                           #
#   o Version for Shiny app development                                 #
# Nov-14-2013							                                              #
#   o Version 0.1 for R package development			                        #
# May-12-2014                                                           #
#   o Version 0.2 - start to incorporate high dimensional case          #
# Sep-23-2014                                                           #
#   o Version 0.3 - different penalties for prognostic coeffs and       #
#                   interaction coeffs, no p-values or CIs for model    #
#                   and regression coeffs in the case of var selection, #
#                   implement 'print' method for eval.pact.cv           #
# Jan-1-2015  and completed in Aug-21-2015                                                          #
#   o Version 0.4 - corrected bug for univar case with nsig=1, show     #
#                   x-axis labels in eval.pact.cv(),                    #
#                   allow for user specified clinical covariates to     #
#                   always be there in the model                        #
# Mar-11-2016                                                           #
#   o version 0.5.0 - Version released to CRAN                          #
#########################################################################


#' @title Predictive Analysis of Clinical Trials
#'
#' @description The \code{pact} package implements a prediction-based approach to the 
#' analysis of data from randomized clinical trials (RCT). Based on clincial response and 
#' covariate data from a RCT comparing a new experimental treatment E versus a control C, 
#' the purpose behind the functions in \code{pact} is to develop and internally validate
#' a model that can identify subjects likely to benefit from E rather than C. Currently, 
#' 'survival' and 'binary' response types are permitted. 
#' 
#'@details
#' Package: pact
#' \cr Type: Package
#' \cr Version: 0.5.0
#' \cr Date: 2016-04-14
#' \cr Author: Dr. Jyothi Subramanian and Dr. Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' \cr License: GPL-3
#' \cr\cr \code{pact.fit} fits a predictive model to data from RCT. Currently, 'survival' and 
#' 'binary' response types are supported. Analysis of high dimensional covariate data is supported.
#' If known and available, a limited number of prognostic covariates can also be specified and 
#' fixed to remain in the predictive model. An object of class 'pact' is returned. 
#' \code{print}, \code{summary} and \code{predict} methods are available for objects of 
#' class 'pact'.
#' Additionally, the function \code{pact.cv} takes as an input the object returned by \code{pact.fit} 
#' and computes predictive scores for each subject through k-fold cross-validation. 
#' Evaluations of the cross-validated predictions are performed by the function \code{eval.pact.cv}. 
#' \cr\cr  Finally, the function \code{overall.analysis}
#' also takes an object of class 'pact' as input and computes some summary statistics 
#' for the comparison of treatments E and C.
#' @docType package
#' @name pact
#' @examples
#' ### Survival response
#' set.seed(10)    
#' data(prostateCancer)     
#' Y <- prostateCancer[,3:4]
#' Xf <- prostateCancer[,7:8]  ## Prognostic covariates fixed to always be in the model
#' Xv <- prostateCancer[,c(5:6,9)]
#' Treatment <- as.factor(prostateCancer[,2])
#' p <- pact.fit(Y=Y,Xf=Xf,Xv=Xv,Treatment=Treatment,family="cox",varSelect="univar")
#' print(p)
#' overall.analysis(p)
#' cv <- pact.cv(p, nfold=5)
#' eval.pact.cv(cv, method="continuous", plot.score=TRUE, perm.test=FALSE, nperm=100)
#' 
#' ### Binary response
#' set.seed(10)
#' data(EORTC10994)
#' Y <- as.factor(EORTC10994[,4])
#' ## No prognostic covariates (Xf) specified
#' Xv <- EORTC10994[,c(2,5:7)]
#' Treatment <- as.factor(EORTC10994[,3])
#' p <- pact.fit(Y=Y,Xv=Xv,Treatment=Treatment,family="binomial",varSelect="none")
#' print(p)
#' overall.analysis(p)
#' cv <- pact.cv(p, nfold=5)
#' eval.pact.cv(cv, method="discrete", g=log(1), perm.test=FALSE, nperm=100)
#' 
#' ### High dimensional data, survival response
#' \dontrun{
#' set.seed(10)    
#' data(GSE10846)     
#' Y <- GSE10846[,1:2]
#' Xv <- GSE10846[,-c(1:3)]
#' Treatment <- as.factor(GSE10846[,3])
#' p <- pact.fit(Y=Y,Xv=Xv,Treatment=Treatment,family="cox",varSelect="lasso",penalty.scaling=2)
#' print(p)
#' overall.analysis(p)
#' cv <- pact.cv(p, nfold=5)
#' eval.pact.cv(cv, method="continuous", plot.score=TRUE, perm.test=FALSE)
#' }
#'
NULL

#' @title Fits a predictive model to the full dataset
#'
#' @description
#' \code{pact.fit} Fits a predictive model using data on all subjects. Currently supports  
#' Cox PH and logistic regression models for 'survival' and 'binary' response types 
#' respectively.
#' 
#' @details
#' A Cox proportional hazards (PH) or a logistic regression model is 
#' developed for data with survival and binary response respectively. Data from subjects 
#' in both 'experimental' (E) and 'control' (C) groups from a RCT is used for model
#' development. Main effect of treatment, main effect of prognostic covariates, 
#' main effects and treatment by covariate interaction terms of candidate predictive  
#' covariates are considered in the model. Methods for variable selection can be optionally   
#' specified by the user for candidate predictive covariates (useful for high-dimensional  
#' covariates). Current options for variable selection include "univar" and "lasso". 
#' In the case of "univar", the number of predictive covariates (nsig) to be included in the
#' model is specified by the user. A univariate selection procedure is applied to identify 
#' covariates that have the lowest treatment*covariate interaction p-values. The 
#' predictive model is then developed using the main effect of treatment, main effects 
#' of prognostic covariates, main effects of the nsig predictive covariates and treatment 
#' by covariate interaction terms for nsig predictive covariates.
#' 
#' In the case of "lasso", an internal cross-validation loop
#' is used to find the penalty value that minimizes the cross-validated error. The user can 
#' choose either the value of the penalty 'lambda' as the penalty that minimizes the 
#' cross-validated  error ("lambda.min") or the largest penalty for which the cross-validated  
#' error is within 1 standard error of the minimum ("lambda.1se"). Also, in the case of "lasso",
#' differential shrinkage can be specified for main effect and interaction effect predictive 
#' coefficients by specifying a value for the ratio of shrinkage for main coefficients to 
#' shrinkage for interaction coefficients. Internally, 'lambda' is scaled using this ratio to 
#' allow for the differential shrinkage of main and interaction coefficients.The penalty 
#' factors affect only variables in Xv and not Xf.
#' 
#' @param Y Response variable. For \code{family='binomial'}, Y should be a factor with two 
#' levels. For \code{family='cox'}, Y should be a two-column matrix with columns named 'time' 
#' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' 
#' indicating right censored. 
#' @param Xf An optional dataframe of the prognostic covariates that are not subject to variable
#' selection and always fixed to remain in the model. Default is NULL (no variable).
#' @param Xv The main dataframe of covariates that are to be used for predictive model development. 
#' Variable selection options affect only the variables in Xv. Xv cannot be NULL.
#' @param Treatment The treatment assignment indicator. A factor with two levels.
#' '0' indicating control (C) and '1' indicating experimental (E) treatment.
#' @param family Type of the response variable. See above. Possible values are 'binomial'
#' or 'cox'.
#' @param varSelect The variable selection method. Possible values are "none","univar" or "lasso".
#' @param nsig The number of covariates to use in the model for varSelect="univar". 
#' Defaults to 3 if the number of candidate covariates is less than 10, else defaults to 10.
#' @param cvfolds.varSelect The number of folds in the internal cross-validation loop for
#' variable selection with varSelect="lasso". Default is 5.
#' @param which.lambda Used with variable selection with varSelect="lasso". Defaults to "min" 
#' if the number of candidate covariates is less than 10, else defaults to "1se". See Details.
#' @param penalty.scaling Ratio of shrinkage applied for main coefficients to shrinkage applied 
#' for interaction coefficients. Used with varSelect="lasso". Default is 0.5. See Details.
#' 
#' @return An object of class 'pact' which is a list with the following components:
#' @return \item{reg}{The fitted regression model}
#' @return \item{family}{Type of the response variable}
#' @return \item{Y}{The response variable used}
#' @return \item{Xf}{The dataframe of prognostic covariates}
#' @return \item{Xv}{The dataframe of candidate predictive variables}
#' @return \item{Treatment}{The treatment assignment indicator used}
#' @return \item{nCovarf}{The number of variables in Xf} 
#' @return \item{nCovarv}{The number of variables in Xv} 
#' @return \item{varSelect}{The variable selection method used}
#' @return \item{nsig, cvfolds.varSelect, which.lambda, penalty.scaling}{The variable selection parameters used}
#' @return \item{call}{The call that produced the return object}
#' 
#' @keywords pact, pact.fit
#' @export
#' @import graphics stats survival glmnet
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' Xf <- prostateCancer[,7:8]
#' Xv <- prostateCancer[,c(5:6,9)]
#' Treatment <- as.factor(prostateCancer[,2])
#' pact.fit(Y=Y, Xf=Xf, Xv=Xv, Treatment=Treatment, family="cox",varSelect="univar")

# Fits a predictive model to the data set and returns an object of class 'pact', 
# which is actually a list 

pact.fit <- function(Y, Xf=NULL, Xv, Treatment, family=c("binomial", "cox"),
                     varSelect=c("none","univar","lasso"),  
                     nsig=ifelse(varSelect == "univar", ifelse(nCovarv < 10,3,10), NA),
                     cvfolds.varSelect=ifelse(varSelect == "lasso", 5, NA), 
                     which.lambda=ifelse(varSelect == "lasso", ifelse(nCovarv < 10,"min","1se"), NA),
                     penalty.scaling=ifelse(varSelect == "lasso", 0.5, NA)) {
  family <- match.arg(family)
  varSelect <- match.arg(varSelect)
  this.call <- match.call()
  
  dimy <- dim(as.matrix(Y))
  nobs <- dimy[1]
  if (nobs == 1) stop("You really want to develop a model with data on just 1 subject?!")
  
  if(is.null(Xv))
    stop("Xv cannot be NULL")
  if (!inherits(Xv, "data.frame")) 
    stop("Xv should be a data frame")

  dimxv <- dim(Xv)
  nrowxv <- ifelse(is.null(dimxv), length(Xv), dimxv[1])
  nCovarv <- ifelse(is.null(dimxv), 1, dimxv[2])
  if (nrowxv != nobs)
    stop(paste("number of subjects in Y (", nobs, ") not equal to the number of subjects in Xv (", 
               nrowxv, ")", sep = ""))

  if (!is.null(Xf)){
    if (!inherits(Xf, "data.frame")) 
      stop("Xf should be a data frame")
    
    dimxf <- dim(Xf)
    nrowxf <- ifelse(is.null(dimxf), length(Xf), dimxf[1])
    nCovarf <- ifelse(is.null(dimxf), 1, dimxf[2])
    if (nrowxf != nobs)
    stop(paste("number of subjects in Y (", nobs, ") not equal to the number of subjects in Xf (", 
               nrowxf, ")", sep = ""))
  } else {
    nCovarf <- 0
  }
  
  Treatment <- as.factor(Treatment)  
  if (!all(match(levels(Treatment), c(0,1),0))) 
    stop("'Treatment' should be a factor with two levels '0' and '1'. See Usage.")
  if (length(Treatment) != nobs)
    stop(paste("number of subjects with Treatment info (", length(Treatment), ") 
               not equal to the number of subjects in Y (", nobs, ")", sep = ""))
  
  res <- switch(family, 
                  cox = .pact.fit.survival(Y, Xf, Xv, Treatment, nCovarf, nCovarv, varSelect, nsig, 
                                           cvfolds.varSelect, penalty.scaling),
                  binomial = .pact.fit.binary(Y, Xf, Xv, Treatment, nCovarf, nCovarv, varSelect, nsig, 
                                              cvfolds.varSelect, penalty.scaling))
  res$Y <- Y
  res$Xf <- Xf
  res$Xv <- Xv
  res$Treatment <- Treatment
  res$nCovarf <- nCovarf
  res$nCovarv <- nCovarv
  res$varSelect <- varSelect
  res$nsig <- nsig
  res$cvfolds.varSelect <- cvfolds.varSelect
  res$which.lambda <- which.lambda
  res$penalty.scaling <- penalty.scaling
  res$call <- this.call
  class(res) <- "pact"
  res
}
 
### Internal function. 'pact.fit' for the "cox" family.

.pact.fit.survival <- function(Y, Xf, Xv, Treatment, nCovarf, nCovarv, varSelect, nsig, cvfolds.varSelect, 
                               penalty.scaling) {
  if (!all(match(c("time", "status"), dimnames(Y)[[2]], 0))) 
    stop("Cox model requires a matrix with columns 'time' (>0) and 'status'  (binary) 
         as a response")
  if (any(Y[, "time"] <= 0)) 
    stop("non-positve event times encountered; not permitted for Cox family")
  
  nEff <- sum(Y[,"status"]) ## Effective sample size is the no. of events (addded check Aug-2015)
  if (nCovarf > nEff/20)
    warning("Seems there are too many fixed prognostic covariates. Proceed with caution!")
  
  SurvObj <- Surv(Y[,"time"], Y[,"status"])
  T <- Treatment

  if (varSelect == "none") {  ## No variable selection 
    cat("No variable selection: All variables present in Xf and Xv used in the model...\n")
    tag <- 0
    if (is.null(Xf))
      tryCatch(CoxReg <- coxph(SurvObj ~ T + . + T:., Xv, na.action="na.exclude"), 
             error=function(err) tag <<- 1)
    else {
      mf1 <- model.frame(~ T + . + T:., Xv, na.action="na.exclude")
      mf2 <- model.frame(~ ., Xf, na.action="na.exclude")
      mm1 <- model.matrix(~ T + . + T:., mf1)
      mm2 <- model.matrix(~ ., mf2)
      mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms
      tryCatch(CoxReg <- coxph(SurvObj ~ ., mm, na.action="na.exclude"), 
               error=function(err) tag <<- 1)
    }
    if (tag == 1) stop("Predictive model cannot be developed using the given dataset.")
  }
  
  if (varSelect == "univar") {  ## univariate variable selection
    cat("Variable selection: Univariate method...\n")
    
    if (nCovarv < nsig)  ## check added 20Aug2015
      stop(paste("number of candidate covariates in Xv is", nCovarv, " and nsig is", nsig))
    
    pval <- .varSelect.uni(Y, Xf, Xv, Treatment, nCovarf, nCovarv, family="cox")
    Covar <- Xv[,order(pval)[1:nsig],drop=FALSE]   ## added 'drop=FALSE'(Jyothi, Jan 2015)
    tag <- 0
    if (is.null(Xf))
      tryCatch(CoxReg <- coxph(SurvObj ~ T + . + T:., Covar, na.action="na.exclude"), 
             error=function(err) tag <<- 1)
    else {
      mf1 <- model.frame(~ T + . + T:., Covar, na.action="na.exclude")
      mf2 <- model.frame(~ ., Xf, na.action="na.exclude")
      mm1 <- model.matrix(~ T + . + T:., mf1)
      mm2 <- model.matrix(~ ., mf2)
      mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms 
      tryCatch(CoxReg <- coxph(SurvObj ~ ., mm, na.action="na.exclude"), 
               error=function(err) tag <<- 1)
    }
    if (tag == 1) stop("Predictive model cannot be developed with this variable selection. Try 'lasso'.")
  } 
  
  if (varSelect == "lasso") {  ## lasso selection
    cat("Variable selection: lasso...\n")
    if (is.null(Xf)) {
      mf <- model.frame(~ T+.+T:., Xv)
      mm <- model.matrix(~ T+.+T:., mf)
      mm <- mm[,-1] ## B'coz we need only the X matrix 
      penalty <- c(0, rep(penalty.scaling,nCovarv), rep(1,nCovarv)) 
    } else {
      mf1 <- model.frame(~ T + . + T:., Xv, na.action="na.exclude")
      mf2 <- model.frame(~ ., Xf, na.action="na.exclude")
      mm1 <- model.matrix(~ T + . + T:., mf1)
      mm2 <- model.matrix(~ ., mf2)
      mm <- cbind(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms
      penalty <- c(0, rep(penalty.scaling,nCovarv), rep(1,nCovarv), rep(0,nCovarf)) 
    }
    ## 'Treatment' is always in the model. Ratio of shrinkage for main coeffs to interaction 
    ## coeffs is specified by the user defined parameter, penalty.scaling (Jyothi, sep 2014)
    ## no penalty for clinical covariates, Xf (Jyothi, Jan 2015)
    
    CoxReg <- glmnet(mm, SurvObj, family="cox", alpha=1, penalty.factor=penalty)
    tag <- 0
    tryCatch(cv <- cv.glmnet(mm, SurvObj, family="cox", alpha=1, penalty.factor=penalty, 
                       nfolds=cvfolds.varSelect), error=function(err) tag <<- 1)
    if (tag == 1) stop("Predictive model cannot be developed with variable selection 'lasso'. Try 'univariate'.")
    CoxReg$lambda.min <- cv$lambda.min
    CoxReg$lambda.1se <- cv$lambda.1se 
  }
  res <- list(reg=CoxReg, family="cox")
  res
  } 

### Internal function. 'pact.fit' for the "binomial" family. 

.pact.fit.binary <- function(Y, Xf, Xv, Treatment, nCovarf, nCovarv, varSelect, nsig, cvfolds.varSelect,
                             penalty.scaling) {  
  dimy = dim(as.matrix(Y))
  nc <- dimy[2]
  
  if (nc == 1) {
    Response = as.factor(Y)
    ntab = table(Y)
    classnames = names(ntab)
  }
  else {
    stop("For binomial family, Y must be a vector with two levels, other dimensions for Y currently not supported")
  }
  
  nEff <- min(ntab)  ## Effective sample size is the no. of obs in the smaller catergory (addded check Aug-2015)
  if (nCovarf > nEff/20)
    warning("Seems there are too many fixed prognostic covariates. Proceed with caution!")
  
  T <- Treatment
  
  if (varSelect == "none") {  ## No variable selection 
    cat("No variable selection: All variables in X used in the model...\n")
    tag <- 0
    if (is.null(Xf))
      tryCatch(LogReg <- glm(Response ~ T + . + T:., Xv, family = binomial(link = "logit"), na.action="na.exclude")
             , error=function(err) tag <<- 1)
    else {
      mf1 <- model.frame(~ T + . + T:., Xv, na.action="na.exclude")
      mf2 <- model.frame(~ ., Xf, na.action="na.exclude")
      mm1 <- model.matrix(~ T + . + T:., mf1)
      mm2 <- model.matrix(~ ., mf2)
      mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercepts
      
      tryCatch(LogReg <- glm(Response ~ ., mm, family = binomial(link = "logit"), na.action="na.exclude")
             , error=function(err) tag <<- 1)
    }
    if (tag == 1) stop("Predictive model cannot be developed using the given dataset.")
  }
  
  if (varSelect == "univar") {  ## univariate variable selection
    cat("Variable selection: Univariate method...\n")
    
    if (nCovarv < nsig)  ## check added 20Aug2015
      stop(paste("number of candidate covariates in Xv is", nCovarv, " and nsig is", nsig))
    
    pval <- .varSelect.uni(Y, Xf, Xv, Treatment, nCovarf, nCovarv, family="binomial")
    Covar <- Xv[,order(pval)[1:nsig],drop=FALSE]   ## added 'drop=FALSE'(Jyothi, Jan 2015)
    tag <- 0
    if (is.null(Xf))
      tryCatch(LogReg <- glm(Response ~ T + . + T:., Covar, family = binomial(link = "logit"), na.action="na.exclude"), 
             error=function(err) tag <<- 1)
    else{
      mf1 <- model.frame(~ T + . + T:., Covar, na.action="na.exclude")
      mf2 <- model.frame(~ ., Xf, na.action="na.exclude")
      mm1 <- model.matrix(~ T + . + T:., mf1)
      mm2 <- model.matrix(~ ., mf2)
      mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## intercept will be fitted by model
      tryCatch(LogReg <- glm(Response ~ ., mm, family = binomial(link = "logit"), na.action="na.exclude"), 
               error=function(err) tag <<- 1)
    }
    if (tag == 1) stop("Predictive model cannot be developed with this variable selection. Try 'lasso'.")
  } 
  
  if (varSelect == "lasso") {  ## lasso selection
    cat("Variable selection: lasso...\n")
    if (is.null(Xf)) {
      mf <- model.frame(~ T+.+T:., Xv)
      mm <- model.matrix(~ T+.+T:., mf)
      mm <- mm[,-1] ## B'coz we need only the X matrix 
      penalty <- c(0, rep(penalty.scaling,nCovarv), rep(1,nCovarv)) 
    } else {
      mf1 <- model.frame(~ T + . + T:., Xv, na.action="na.exclude")
      mf2 <- model.frame(~ ., Xf, na.action="na.exclude")
      mm1 <- model.matrix(~ T + . + T:., mf1)
      mm2 <- model.matrix(~ ., mf2)
      mm <- cbind(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercepts
      penalty <- c(0, rep(penalty.scaling,nCovarv), rep(1,nCovarv), rep(0,nCovarf)) 
    }
    ## 'Treatment' is always in the model. Ratio of shrinkage for main coeffs to interaction 
    ## coeffs is specified by the user defined parameter, penalty.scaling (Jyothi, sep 2014)
    ## no penalty for clinical covariates, Xf (Jyothi, Jan 2015)
    
    LogReg <- glmnet(mm, Response, family="binomial", alpha=1, penalty.factor=penalty)
    tag <- 0
    tryCatch(cv <- cv.glmnet(mm, Response, family="binomial", alpha=1, penalty.factor=penalty, 
                                     nfolds=cvfolds.varSelect), error=function(err) tag <<- 1)
    if (tag == 1) stop("Predictive model cannot be developed with variable selection 'lasso'. Try 'univariate'.")
    LogReg$lambda.min <- cv$lambda.min
    LogReg$lambda.1se <- cv$lambda.1se 
  }
  res <- list(reg=LogReg, family="binomial")
  res
}

#' @title Print an object of class 'pact' 
#'
#' @description
#' print method for objects of class 'pact'
#' 
#' @details
#' The call that produced the object is printed, followed by the classification 
#' function from \code{pact.fit} for calculating the predictive scores for new 
#' subjects
#' 
#' @method print pact
#' 
#' @param x The object returned from \code{'pact.fit'}
#' @param digits significant digits in the print 
#' @param ... Additional print arguments
#' 
#' @return The classification function is printed
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' Xf <- prostateCancer[,7:8]
#' Xv <- prostateCancer[,c(5:6,9)]
#' Treatment <- as.factor(prostateCancer[,2])
#' p <- pact.fit(Y=Y, Xf=Xf, Xv=Xv, Treatment=Treatment, family="cox", varSelect="lasso")
#' print(p)

print.pact <- function(x, digits = max(3, getOption("digits") - 3), ...) {	### x = object of class 'pact'
  reg <- x$reg
  varSelect <- x$varSelect
  nCovarf <- x$nCovarf
  nCovarv <- x$nCovarv
  cat("\nCall: ", deparse(x$call), "\n\n")
  cat("\nfamily: ", x$family, "\n\n")
  
  if (varSelect != "lasso") {
    if (x$family == "cox") {
      nInteractions <- (nrow(summary(reg)$coef) - 1 - nCovarf)/2
      ClassFn <- paste("(", round(summary(reg)$coef[1,1], digits), ")")
      for (i in 2:(nInteractions + 1)) {
        ClassFn <- paste(ClassFn, "\n", "  +  (", round(summary(reg)$coef[nInteractions+i, 1], digits), " )*", row.names(summary(reg)$coef)[i], sep = "")
      }
      cat(paste("\nClassification function for classifying future subjects:\n", "f = ", ClassFn, "\n\n"), sep = "")
    } else if (x$family == "binomial") { ### For binary response (intercept also exists in the model output for binary response)
      nInteractions <- (nrow(summary(reg)$coef) - 2 - nCovarf)/2
      ClassFn <- paste("(", round(summary(reg)$coef[1,1], digits), ")")
      ClassFn <- paste(ClassFn, "\n", "  +  (", round(summary(reg)$coef[2,1], digits), " )")
      for (i in 3:(nInteractions + 2)) {
        ClassFn <- paste(ClassFn, "\n", "  +  (", round(summary(reg)$coef[nInteractions+i, 1], digits), " )*", row.names(summary(reg)$coef)[i], sep = "")
      }
      cat(paste("\nClassification function for classifying future subjects:\n", 
                "f = ", ClassFn, "\n\n"), sep = "")
    }
  } else {  ### if varSelect == "lasso"
    s <- ifelse(x$which.lambda == "min", reg$lambda.min, reg$lambda.1se)
    coeffs <- coef(reg, s = s)
    if (x$family == "cox") {
      nInteractions <- (nrow(coeffs) - 1 - nCovarf)/2
      ClassFn <- paste("(", round(coeffs[1], digits), ")")
      for (i in 2:(nInteractions + 1)) {
        if (coeffs[nInteractions + i] != 0) 
          ClassFn <- paste(ClassFn, "\n", "  +  (", round(coeffs[nInteractions+i], digits), " )*", row.names(coeffs)[i], sep = "")
      }
      cat(paste("\nClassification function for classifying future subjects:\n", "f = ", ClassFn, "\n\n"), sep = "")
    }
    else if (x$family == "binomial") {
      nInteractions <- (nrow(coeffs) - 2 - nCovarf)/2
      ClassFn <- paste("(", round(coeffs[1], digits), ")")
      ClassFn <- paste(ClassFn, "\n", "  +  (", round(coeffs[2], digits), " )")
      for (i in 3:(nInteractions + 2)) {
        if (coeffs[nInteractions + i] != 0) 
          ClassFn <- paste(ClassFn, "\n", "  +  (", round(coeffs[nInteractions+i], digits), " )*", row.names(coeffs)[i], sep = "")
      }
      cat(paste("\nClassification function for classifying future subjects:\n", "f = ", ClassFn, "\n\n"), sep = "")
    }
  }  
}

#' @title Summarize a predictive model fit
#'
#' @description
#' summary method for objects of class 'pact'
#' 
#' @details
#' Returns all coefficient estimates from the regression model of the 'pact' object 
#' 
#' @method summary pact
#' 
#' @param object The object returned from \code{'pact.fit'}
#' @param ... Additional arguments for 'summary'
#' 
#' @return All the coefficient estimates from the regression model fitted by \code{pact.fit}
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' Xf <- prostateCancer[,7:8]
#' Xv <- prostateCancer[,c(5:6,9)]
#' Treatment <- as.factor(prostateCancer[,2])
#' p <- pact.fit(Y=Y,Xf=Xf,Xv=Xv,Treatment=Treatment,family="cox",varSelect="none")
#' summary(p)
#' 
summary.pact <- function(object, ...) {	### object = object of class 'pact'
  reg <- object$reg
  family <- object$family
  varSelect <- object$varSelect
  
  if (varSelect != "lasso") {
     coeffs <- summary(reg)$coeff[,1]
  } else {  ## varSelect == "lasso"
    a <- ifelse(object$which.lambda == "min", reg$lambda.min, reg$lambda.1se)
    coeffs <- coef(reg, s=a)
    coeffs <- coeffs[which(coeffs[,1] != 0),,drop=FALSE]
  }
  coeffs
}

#' @title Predictions from a predictive model fit
#'
#' @description
#' Predicts the scores for new subjects from a previously developed object of class 'pact'
#' 
#' @details
#' Returns the scores for new subjects from an object of class 'pact', given their covariate values and 
#' treatment assignment. 
#' 
#' @method predict pact
#' 
#' @param object The object returned from \code{pact.fit} 
#' @param newxv The dataframe Xv of covariates for the new subjects for whom predictions 
#' are to be made
#' @param ... Other arguments to 'predict'
#' 
#' @return A numeric vector containing the predicted scores for the new subjects 
#' from the fitted model is returned
#'
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' ### Survival response
#' data(prostateCancer)
#' Y <- prostateCancer[1:400,3:4]
#' Xf <- prostateCancer[1:400,7:8]
#' Xv <- prostateCancer[1:400,c(5:6,9)]
#' Treatment <- as.factor(prostateCancer[1:400,2])
#' p <- pact.fit(Y=Y, Xf=Xf, Xv=Xv, Treatment=Treatment, family="cox", varSelect="univar")
#' 
#' newxv <- prostateCancer[401:410,c(5:6,9)]
#' predict(p, newxv)
#' 
#' ### Binary response
#' data(EORTC10994)
#' Y <- as.factor(EORTC10994[1:120,4])
#' Xv <- EORTC10994[1:120,c(2,5:7)]
#' Treatment <- as.factor(EORTC10994[1:120,3])
#' p <- pact.fit(Y=Y,Xv=Xv,Treatment=Treatment,family="binomial",varSelect="none")
#' 
#' newxv <- EORTC10994[121:125,c(2,5:7)]
#' predict(p, newxv)

predict.pact <- function(object, newxv, ...) {  ### object = object of class 'pact'
  if (missing(newxv)) {
    stop("You need to supply a value for 'newxv'")
  }
  if (!inherits(newxv, "data.frame")) 
    stop("newxv should be a data frame")
  
  dimxv <- dim(newxv)
  nobs <- dimxv[1]
  nCovarv <- dimxv[2]
  
  ## Can add a check here for the presence of the same variables as in developed model...(26May14)
  ## No need..check shd be there in the original predict method...(29May2014)
  
  reg <- object$reg
  family <- object$family
  varSelect <- object$varSelect
  Xf <- object$Xf
  
  T <- factor(rep(1, nobs),levels=c(0,1))
  
  if (varSelect != "lasso") {
    type <- ifelse(family == "cox","lp","link")
    if(is.null(Xf)) {
      LinPred <- predict(reg, cbind(T,newxv), type=type)      
        #### Now predicting for inverted treatment assignment
      T <- factor(rep(0, nobs), levels=c(0,1))
      InvLinPred <- predict(reg, cbind(T,newxv), type=type)
    }
    else {
      newxf <- Xf[rep(1,each=nobs),,drop=FALSE]
      mf1 <- model.frame(~ T + . + T:., newxv, na.action="na.exclude")
      mf2 <- model.frame(~ ., newxf, na.action="na.exclude")
      mm1 <- model.matrix(~ T + . + T:., mf1)
      mm2 <- model.matrix(~ ., mf2)
      mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept
      LinPred <- predict(reg, mm, type=type)
        #### Now predicting for inverted treatment assignment
      T <- factor(rep(0, nobs), levels=c(0,1))
      mf1 <- model.frame(~ T + . + T:., newxv, na.action="na.exclude")
      mf2 <- model.frame(~ ., newxf, na.action="na.exclude")
      mm1 <- model.matrix(~ T + . + T:., mf1)
      mm2 <- model.matrix(~ ., mf2)
      mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept
      InvLinPred <- predict(reg, mm, type=type)
    }
      #### (log HR/odds for T=1 - log HR/odds for T=0)
    PredScore <- LinPred - InvLinPred
  }
  else {  ### varSelect == "lasso"
    s <- ifelse(object$which.lambda == "min", reg$lambda.min, reg$lambda.1se)
    if(is.null(Xf)){
      mf <- model.frame(~ T+.+T:., newxv)
      mm <- model.matrix(~ T+.+T:., mf)
      mm <- mm[,-1]
      LinPred <- predict(reg, mm, s=s, type="link")
    #### Now predicting for inverted treatment assignment
      T <- factor(rep(0, nobs), levels=c(0,1))
      mf <- model.frame(~ T+.+T:., newxv)
      mm <- model.matrix(~T+.+T:., mf)
      mm <- mm[,-1]
      InvLinPred <- predict(reg, mm, s=s, type="link")
    } else {
      newxf <- Xf[rep(1,each=nobs),,drop=FALSE]
      mf1 <- model.frame(~ T + . + T:., newxv, na.action="na.exclude")
      mf2 <- model.frame(~ ., newxf, na.action="na.exclude")
      mm1 <- model.matrix(~ T + . + T:., mf1)
      mm2 <- model.matrix(~ ., mf2)
      mm <- cbind(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept
      LinPred <- predict(reg, mm, s=s, type="link")
     #### Now predicting for inverted treatment assignment
      T <- factor(rep(0, nobs), levels=c(0,1))
      mf1 <- model.frame(~ T + . + T:., newxv, na.action="na.exclude")
      mf2 <- model.frame(~ ., newxf, na.action="na.exclude")
      mm1 <- model.matrix(~ T + . + T:., mf1)
      mm2 <- model.matrix(~ ., mf2)
      mm <- cbind(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept
      InvLinPred <- predict(reg, mm, s=s, type="link")
    }
      #### (log HR/odds for T=1 - log HR/odds for T=0)
    PredScore <- LinPred - InvLinPred
  }
  PredScore
}


### End of main 'pact' fucntions
### Other functions

### Evaluation of an object of class 'pact' using cross-validation

#' @title Split a dataset into k parts for k-fold cross-validation
#'
#' @description
#' Split a dataset into k parts for k-fold cross-validation. This function is used in 
#' \code{pact.cv} to create the splits for cross-validation
#' 
#' @param n The sample size 
#' @param k The number of folds. k=n would mean a leave-one-out cross-validation
#'  
#' @return A integer vector of same length as n. Each observation is an integer taking a value 
#' between 1 to k denoting the fold it belongs to.
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' KfoldCV(15,3)
#' KfoldCV(15,15)
 
KfoldCV <- function(n,k) { ### Partition a dataset to k parts for k-fold CV
  f <- ceiling(n/k)
  s <- sample(rep(1:k, f), n)
  s
}


### Internal function. Univariate variable selection for the high dim case
.varSelect.uni <- function(Y, Xf, Xv, Treatment, nCovarf, nCovarv, family = c("binomial", "cox")) {
  Covar <- Xv
  T <- Treatment
  p.uni <- NULL
  
  if (family == "cox") {
    SurvObj <- Surv(Y[,"time"], Y[,"status"])
    pos <- 3  ## T+var+T*var+Clincovar (modified Jan2015)
    for (i in 1:nCovarv) {
      tag <- 0
      if (is.null(Xf)) {
        tryCatch(t <- coxph(SurvObj ~ T + . + T:., data = data.frame(Covar[, i])), 
                error=function(err) tag <<- 1)
      } else {
        mf1 <- model.frame(~ T + . + T:., data.frame(Covar[, i]), na.action="na.exclude")
        mf2 <- model.frame(~ ., Xf, na.action="na.exclude")
        mm1 <- model.matrix(~ T + . + T:., mf1)
        mm2 <- model.matrix(~ ., mf2)
        mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms 
        tryCatch(t <- coxph(SurvObj ~ ., mm, na.action="na.exclude"), 
                 error=function(err) tag <<- 1)
      }
      ifelse((tag == 0), p.uni[i] <- summary(t)$coef[pos,5], p.uni[i] <- 1)
    }
    } else {	## family = "binomial"
    Y = as.factor(Y)
    pos <- 4	## Intercept+T+var+T*var+ClinCovar(if any) (modified Jan2015)
    for (i in 1:nCovarv) {
      tag <- 0
      if (is.null(Xf)) {
        tryCatch(t <- glm(Y ~ T + . + T:. ,data = data.frame(Covar[, i]), family = binomial(link = "logit")),
               error=function(err) tag <<- 1)
      } else {
        mf1 <- model.frame(~ T + . + T:., data.frame(Covar[, i]), na.action="na.exclude")
        mf2 <- model.frame(~ ., Xf, na.action="na.exclude")
        mm1 <- model.matrix(~ T + . + T:., mf1)
        mm2 <- model.matrix(~ ., mf2)
        mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept
        tryCatch(t <- glm(Response ~ ., mm, family = binomial(link = "logit"), na.action="na.exclude"), 
                 error=function(err) tag <<- 1)
       }
      ifelse((tag == 0), p.uni[i] <- summary(t)$coef[pos,4], p.uni[i] <- 1)
    }
  }
  p.uni
}


#' @title Cross-validation for pact
#'
#' @description
#' Predictive scores using k-fold cross-validation for the model developed in \code{pact.fit}
#' 
#' @details
#' Obtain cross-validated predictive scores for the model developed in \code{pact.fit}.
#' In each fold of the cross-validation, a model is developed from the observations in the 
#' training set using the same variable selection parameters as that used for the model  
#' developed in \code{pact.fit}. The estimated coefficients of the regression model developed 
#' using training set are used to make predictions for the left out observations (test set). 
#' This is repeated for all the folds. Scores are thus obtained for all the subjects in the dataset. 
#' The function \code{\link{eval.pact.cv}} provides various evaluation options for the cross-validated 
#' scores.
#' 
#' @param p An object of class 'pact'
#' @param nfold The number of folds (k) for the k-fold cross-validation. k equal to the sample size 
#' would mean a leave-one-out cross-validation
#' 
#' @return A list with the following components
#' @return \item{PredScore}{The cross-validated scores for each subject (a vector)}
#' @return \item{Y}{The response variable used}
#' @return \item{Xf}{The dataframe of fixed prognostic covariates}
#' @return \item{Xv}{The dataframe of candidate predictive variables}
#' @return \item{Treatment}{The treatment assignment indicator used}
#' @return \item{nCovarf}{The number of variables in Xf} 
#' @return \item{nCovarv}{The number of variables in Xv} 
#' @return \item{family}{Type of the response variable}
#' @return \item{varSelect}{The variable selection method used}
#' @return \item{nsig, cvfolds.varSelect, which.lambda, penalty.scaling}{The variable selection parameters used}
#' @return \item{call}{The call that produced this output}
#' 
#' @keywords pact, pact.cv
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' Xf <- prostateCancer[,7:8]
#' Xv <- prostateCancer[,c(5:6,9)]
#' Treatment <- as.factor(prostateCancer[,2])
#' p <- pact.fit(Y=Y,Xf=Xf,Xv=Xv,Treatment=Treatment,family="cox",varSelect="lasso")
#' cv <- pact.cv(p, nfold=5)

pact.cv <- function(p, nfold) {
  if(!inherits(p,"pact"))
    stop("First argument must be an object of class 'pact'.")
  family <- p$family
  this.call <- match.call()
  varSelect <- p$varSelect
  nsig <- p$nsig
  cvfolds.varSelect <- p$cvfolds.varSelect
  which.lambda <- p$which.lambda
  penalty.scaling <- p$penalty.scaling
  
  Y <- p$Y
  Xf <- p$Xf
  Xv <- p$Xv

  Treatment <- as.factor(p$Treatment)  
  nCovarf <- p$nCovarf
  nCovarv <- p$nCovarv
  
  dimxv <- dim(Xv)
  nobs <- dimxv[1]
  nfold <- as.integer(nfold)
  if ((nfold > nobs) || (nfold < 2))
    stop("nfold should be an integer >= 2 and cannot be greater than the number of observations")
  
  #### K-fold CV ####
  res <- switch(family, 
                cox = .pact.cv.survival(Y, Xf, Xv, Treatment, nCovarf, nCovarv, nfold, varSelect,  
                                        nsig, cvfolds.varSelect, which.lambda, penalty.scaling),
                binomial = .pact.cv.binary(Y, Xf, Xv, Treatment, nCovarf, nCovarv, nfold, varSelect,  
                                           nsig, cvfolds.varSelect, which.lambda, penalty.scaling))
  res$Y <- Y
  res$Xf <- Xf
  res$Xv <- Xv
  res$Treatment <- Treatment
  res$nCovarf <- nCovarf
  res$nCovarv <- nCovarv
  res$nfold <- nfold
  res$varSelect <- varSelect
  res$nsig <- nsig
  res$cvfolds.varSelect <- cvfolds.varSelect
  res$which.lambda <- which.lambda
  res$penalty.scaling <- penalty.scaling
  res$call <- this.call
  class(res) <- "pact.cv"
  res
}

### Internal cross-validation function. For survival response.
.pact.cv.survival <- function(Y, Xf, Xv, Treatment, nCovarf, nCovarv, nfold, varSelect, nsig, 
                              cvfolds.varSelect, which.lambda, penalty.scaling) { 
  SurvObj <- Surv(Y[, "time"],Y[, "status"])
  dimxv <- dim(Xv)
  nobs <- dimxv[1]
  
  #### Generate IDs for K-fold cross-validation
  
  CV <- KfoldCV(nobs, nfold)
  
  #### K-fold CV
  
  Temp.CVest <- lapply(1:nfold,function(i) { 
    Ind.train <- which(CV != i)
    Ind.test <- which(CV == i)
    
    if (varSelect != "lasso") {
      T <- Treatment[Ind.train]
      if (varSelect == "univar") {  ## univariate var selection
        pval <- .varSelect.uni(Y[Ind.train,], Xf[Ind.train,], Xv[Ind.train,], T, nCovarf, nCovarv, family="cox")
        sig.X <- Xv[,order(pval)[1:nsig],drop=FALSE]   ## added 'drop=FALSE'(Jyothi, Jan 2015)
      } else { ## no var selection
        sig.X <- Xv
      }
      tag <- 0
      if (is.null(Xf))
        tryCatch(CoxReg <- coxph(SurvObj[Ind.train] ~ T + . + T:., sig.X[Ind.train,,drop=FALSE], 
                                 na.action="na.exclude"), error=function(err) tag <<- 1)
      else {
        mf1 <- model.frame(~ T + . + T:., sig.X[Ind.train,,drop=FALSE], na.action="na.exclude")
        mf2 <- model.frame(~ ., Xf[Ind.train,,drop=FALSE], na.action="na.exclude")
        mm1 <- model.matrix(~ T + . + T:., mf1)
        mm2 <- model.matrix(~ ., mf2)
        mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms 
        tryCatch(CoxReg <- coxph(SurvObj[Ind.train] ~ ., mm, na.action="na.exclude"), 
                 error=function(err) tag <<- 1)
      }
      if (tag == 0) {
          #### Now predict scores for the test set cases using CoxReg and the treatment inversion steps
        T <- Treatment[Ind.test]
        if (is.null(Xf))
          LinPred.test <- predict(CoxReg, cbind(T,sig.X[Ind.test,,drop=FALSE]), "lp")
        else {
          mf1 <- model.frame(~ T + . + T:., sig.X[Ind.test,,drop=FALSE], na.action="na.exclude")
          mf2 <- model.frame(~ ., Xf[Ind.test,,drop=FALSE], na.action="na.exclude")
          mm1 <- model.matrix(~ T + . + T:., mf1)
          mm2 <- model.matrix(~ ., mf2)
          mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept
          LinPred.test <- predict(CoxReg, mm, "lp")
        }
        T <- as.factor(ifelse((T == 1),0,1))
        if(is.null(Xf))
          InvLinPred.test <- predict(CoxReg, cbind(T,sig.X[Ind.test,,drop=FALSE]), "lp")
        else {
          mf1 <- model.frame(~ T + . + T:., sig.X[Ind.test,,drop=FALSE], na.action="na.exclude")
          mf2 <- model.frame(~ ., Xf[Ind.test,,drop=FALSE], na.action="na.exclude")
          mm1 <- model.matrix(~ T + . + T:., mf1)
          mm2 <- model.matrix(~ ., mf2)
          mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept
          InvLinPred.test <- predict(CoxReg, mm, "lp")
        }
        PredScore.test <- ifelse((T == 0), LinPred.test - InvLinPred.test, InvLinPred.test - LinPred.test)
      } else { #### In case the Cox regression model had failed to converge
        PredScore.test <- NA
      }
    }
    
    if (varSelect == "lasso"){
      T <- Treatment[Ind.train]
      tempX <- Xv[Ind.train,,drop=FALSE]
      
      if (is.null(Xf)) {
        mf <- model.frame(~ T+.+T:., tempX)
        mm <- model.matrix(~ T+.+T:., mf)
        mm <- mm[,-1] ## B'coz we need only the X matrix for glmnet
        penalty <- c(0, rep(penalty.scaling,nCovarv), rep(1,nCovarv)) 
      } else {
        mf1 <- model.frame(~ T + . + T:., tempX, na.action="na.exclude")
        mf2 <- model.frame(~ ., Xf[Ind.train,,drop=FALSE], na.action="na.exclude")
        mm1 <- model.matrix(~ T + . + T:., mf1)
        mm2 <- model.matrix(~ ., mf2)
        mm <- cbind(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms 
        penalty <- c(0, rep(penalty.scaling,nCovarv), rep(1,nCovarv), rep(0,nCovarf)) 
      }
      tag <- 0
      tryCatch(CoxReg <- glmnet(mm, SurvObj[Ind.train], family="cox", alpha=1, penalty.factor=penalty),
               error=function(err) tag <<- 1)
      if (tag == 0) {
        cv <- cv.glmnet(mm, SurvObj[Ind.train], family="cox", alpha=1, penalty.factor=penalty, 
                               nfolds=cvfolds.varSelect)
        s <- ifelse(which.lambda == "min", cv$lambda.min, cv$lambda.1se)
          #### Now predict for the test set cases using CoxReg and by repeating the treatment inversion steps
        tempX <- Xv[Ind.test,,drop=FALSE]
        T <- Treatment[Ind.test]
        if(is.null(Xf)) {
          mf <- model.frame(~ T+.+T:., tempX)
          mm <- model.matrix(~T+.+T:., mf)
          mm <- mm[,-1]
        } else {
          mf1 <- model.frame(~ T + . + T:., tempX, na.action="na.exclude")
          mf2 <- model.frame(~ ., Xf[Ind.test,,drop=FALSE], na.action="na.exclude")
          mm1 <- model.matrix(~ T + . + T:., mf1)
          mm2 <- model.matrix(~ ., mf2)
          mm <- cbind(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms 
        }
        LinPred.test <- predict(CoxReg, mm, s=s, type="link")
        
        T <- as.factor(ifelse((T == 1),0,1))
        if (is.null(Xf)) {
          mf <- model.frame(~ T+.+T:., tempX)
          mm <- model.matrix(~T+.+T:., mf)
          mm <- mm[,-1]
        } else {
          mf1 <- model.frame(~ T + . + T:., tempX, na.action="na.exclude")
          mf2 <- model.frame(~ ., Xf[Ind.test,,drop=FALSE], na.action="na.exclude")
          mm1 <- model.matrix(~ T + . + T:., mf1)
          mm2 <- model.matrix(~ ., mf2)
          mm <- cbind(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms 
        }
        InvLinPred.test <- predict(CoxReg, mm, s=s, type="link")
        PredScore.test <- ifelse((T == 0), LinPred.test - InvLinPred.test, InvLinPred.test - LinPred.test)
      } else { #### In case the Cox regression model had failed to converge
        PredScore.test <- NA
      }
    }
    ifelse((tag == 0), return(PredScore.test), return(NA))
  })
  
  PredScore <- NULL 
  for(i in (1:nfold)) {
    PredScore[CV == i] <- Temp.CVest[[i]]
  }
  out.cv <- list(PredScore=PredScore, family="cox")
  out.cv
}

### Internal cross-validation function. For binary response.

.pact.cv.binary <- function(Y, Xf, Xv, Treatment, nCovarf, nCovarv, nfold, varSelect, nsig, cvfolds.varSelect,
                            which.lambda, penalty.scaling) { 
  Response <- Y
  dimxv <- dim(Xv)
  nobs <- dimxv[1]

  #### Generate IDs for K-fold cross-validation
  
  CV <- KfoldCV(nobs, nfold)
  
  #### K-fold CV
  
  Temp.CVest <- lapply(1:nfold,function(j) { 
    Ind.train <- which(CV != j)
    Ind.test <- which(CV == j)
    
    if (varSelect != "lasso") {
      T <- Treatment[Ind.train]
      if (varSelect == "univar") {  ## univariate var selection
        pval <- .varSelect.uni(Response[Ind.train], Xf[Ind.train,], Xv[Ind.train,], T, nCovarf, nCovarv, family="binomial")
        sig.X <- Xv[,order(pval)[1:nsig],drop=FALSE]   ## added 'drop=FALSE'(Jyothi, Jan 2015)
      } else {  ## no var selection
        sig.X <- Xv
      }
      
      tag <- 0
      if (is.null(Xf))
        tryCatch(LogReg <- glm(Response[Ind.train] ~ T + . + T:., sig.X[Ind.train,,drop=FALSE], 
                             family = binomial(link = "logit"), na.action="na.exclude"), 
                              error=function(err) tag <<- 1)
      else {
        mf1 <- model.frame(~ T + . + T:., sig.X[Ind.train,,drop=FALSE], na.action="na.exclude")
        mf2 <- model.frame(~ ., Xf[Ind.train,,drop=FALSE], na.action="na.exclude")
        mm1 <- model.matrix(~ T + . + T:., mf1)
        mm2 <- model.matrix(~ ., mf2)
        mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms 
        tryCatch(LogReg <- glm(Response[Ind.train] ~ ., mm, family = binomial(link = "logit"), 
                               na.action="na.exclude"), error=function(err) tag <<- 1)
      }
      
      if (tag == 0) {
        #### Now predicting for the test set cases using LogReg and by repeating the treatment inversion steps
        T <- Treatment[Ind.test]
        if (is.null(Xf))
          LinPred.test <- predict(LogReg, cbind(T,sig.X[Ind.test,,drop=FALSE]), type="link")
        else{
          mf1 <- model.frame(~ T + . + T:., sig.X[Ind.test,,drop=FALSE], na.action="na.exclude")
          mf2 <- model.frame(~ ., Xf[Ind.test,,drop=FALSE], na.action="na.exclude")
          mm1 <- model.matrix(~ T + . + T:., mf1)
          mm2 <- model.matrix(~ ., mf2)
          mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept
          LinPred.test <- predict(LogReg, mm, "link")
        }
        T <- as.factor(ifelse((T == 1),0,1))
        if (is.null(Xf))
          InvLinPred.test <- predict(LogReg, cbind(T,sig.X[Ind.test,,drop=FALSE]), type="link")
        else {
          mf1 <- model.frame(~ T + . + T:., sig.X[Ind.test,,drop=FALSE], na.action="na.exclude")
          mf2 <- model.frame(~ ., Xf[Ind.test,,drop=FALSE], na.action="na.exclude")
          mm1 <- model.matrix(~ T + . + T:., mf1)
          mm2 <- model.matrix(~ ., mf2)
          mm <- data.frame(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept
          InvLinPred.test <- predict(LogReg, mm, type="link")
        }
        PredScore.test <- ifelse((T == 0), LinPred.test - InvLinPred.test, InvLinPred.test - LinPred.test)
      } else { #### In case the regression model failed to converge
        PredScore.test <- NA
      }
    }
    
    if (varSelect == "lasso"){
      T <- Treatment[Ind.train]
      tempX <- Xv[Ind.train,,drop=FALSE]
      
      if(is.null(Xf)){
        mf <- model.frame(~ T+.+T:., tempX)
        mm <- model.matrix(~T+.+T:., mf)
        mm <- mm[,-1]
        penalty <- c(0, rep(penalty.scaling,nCovarv), rep(1,nCovarv)) 
      } else {
        mf1 <- model.frame(~ T + . + T:., tempX, na.action="na.exclude")
        mf2 <- model.frame(~ ., Xf[Ind.train,,drop=FALSE], na.action="na.exclude")
        mm1 <- model.matrix(~ T + . + T:., mf1)
        mm2 <- model.matrix(~ ., mf2)
        mm <- cbind(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms 
        penalty <- c(0, rep(penalty.scaling,nCovarv), rep(1,nCovarv), rep(0,nCovarf)) 
      }
      
      tag <- 0
      tryCatch(LogReg <- glmnet(mm, Response[Ind.train], family="binomial", alpha=1, penalty.factor=penalty),
               error=function(err) tag <<- 1)
      if (tag == 0) {
        cv <- cv.glmnet(mm, Response[Ind.train], family="binomial", alpha=1, penalty.factor=penalty, 
                                 nfolds=cvfolds.varSelect)
        s <- ifelse(which.lambda == "min", cv$lambda.min, cv$lambda.1se)
        #### Now predict for the test set cases using CoxReg and by repeating the treatment inversion steps
        tempX <- Xv[Ind.test,,drop=FALSE]
        T <- Treatment[Ind.test]
        if(is.null(Xf)) {
          mf <- model.frame(~ T+.+T:., tempX)
          mm <- model.matrix(~T+.+T:., mf)
          mm <- mm[,-1]
        } else {
          mf1 <- model.frame(~ T + . + T:., tempX, na.action="na.exclude")
          mf2 <- model.frame(~ ., Xf[Ind.test,,drop=FALSE], na.action="na.exclude")
          mm1 <- model.matrix(~ T + . + T:., mf1)
          mm2 <- model.matrix(~ ., mf2)
          mm <- cbind(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms 
        }
        LinPred.test <- predict(LogReg, mm, s=s, type="link")
        
        T <- as.factor(ifelse((T == 1),0,1))
        if(is.null(Xf)) {
          mf <- model.frame(~ T+.+T:., tempX)
          mm <- model.matrix(~T+.+T:., mf)
          mm <- mm[,-1]
        } else {
          mf1 <- model.frame(~ T + . + T:., tempX, na.action="na.exclude")
          mf2 <- model.frame(~ ., Xf[Ind.test,,drop=FALSE], na.action="na.exclude")
          mm1 <- model.matrix(~ T + . + T:., mf1)
          mm2 <- model.matrix(~ ., mf2)
          mm <- cbind(mm1[,-1,drop=FALSE],mm2[,-1,drop=FALSE]) ## remove intercept terms 
        }
        InvLinPred.test <- predict(LogReg, mm, s=s, type="link")
        
        PredScore.test <- ifelse((T == 0), LinPred.test - InvLinPred.test, InvLinPred.test - LinPred.test)
      } else { #### In case the regression model had failed to converge
        PredScore.test <- NA
      }
    }
    ifelse((tag == 0), return(PredScore.test), return(NA))
  })
  
  PredScore <- NULL 
  for(i in (1:nfold)) {
    PredScore[CV == i] <- Temp.CVest[[i]]
  }
  out.cv <- list(PredScore=PredScore, family="binomial")
  out.cv
}


#' @title Evaluation functions for cross-validated predictions
#'
#' @description
#' Methods for the evaluation of the cross-validated predictive scores obtained from \code{\link{pact.cv}}
#' 
#' @details
#' Currently two methods are defined for the evaluation of the scores obtained from \code{pact.cv}. In 
#' \code{method='discrete'} a user specified cut-off score is used to classify the subjects into groups
#' 'benefit' or 'do not benefit' from new treatment. In each of the 'benefit' and 'do not benefit' groups 
#' the actual responses in the control (C) and the experimental (E) groups are compared.   
#' For the 'cox' family, the 'score' for a subject represents the predicted change in the log hazard when 
#' the subject is treated with E as against C (with lower values denoting benefit with E). In the case of the 
#' 'binomial' family, the 'score' represents the predicted change in the log odds of a response when the 
#' subject is treated with E as against C (with higher values denoting benefit with E). 
#' For the 'cox' family, examples of the cut-point \code{g} could be \code{g=log(1)} with score < g
#' meaning benefit with E. Or one could be more stringent and have \code{g} correspond to 
#' a 30\% reduction in hazard (\code{g=log(0.70)}). 
#' For the 'binomial' family, \code{g=log(1.20)} with score > g meaning sensitive to E would mean that
#' subjects predicted to receive at least 20\% increase in odds of response with E are 
#' classified as benefitting from E.  
#' \cr\cr In \code{method='continuous'} no cut-off is applied to the cross-validated scores. 
#' A Cox proportional hazards (PH) regression or a logistic regression 
#' model (respectively for 'survival' and 'binary' response) is then developed that includes 
#' the main effect of treatment, main effect of cross-validated score, and treatment*score interaction.
#' For survival response, this model is used to generate the Kaplan Meier survival curves for each treatment 
#' at the at 20th, 40th, 60th and 80th percentiles of predictive scores (\code{plot.score = TRUE}). 
#' The model is also used to compute the estimated probability of surviving beyond a landmark time 
#' specified in \code{plot.time} as a function of treatment and (cross-validated) score (if 
#' \code{plot.time = NULL}, this plot is not produced). For binary response, the output from evaluation
#' is a plot of the probability of response as a functions of the predictive score and Treatment.
#' \cr\cr If \code{perm.test=TRUE}, permutation based significance tests are performed on appropriate 
#' test statistics and p-values are computed. See 'Value' and the package vignette and for more details 
#' on the permutation tests.
#'
#' @param out.cv The object from \code{pact.cv} 
#' @param method The evaluation method. Currently two options, \code{method='discrete'} or
#' \code{method='continuous'}, are available. See 'Details'.
#' @param g The cut-point for grouping scores into subsets 'benefit' and 'no benefit' from 
#' new treatment. Ignored for \code{method='continuous'}.
#' @param plot.score Used only for plots if \code{method='continuous'} is chosen for survival response. 
#' Logical representing whether survival 
#' curves at specific quantiles of cross-validated scores are to be drawn. See 'Details'.
#' @param plot.time Used only for plots if \code{method='continuous'} is chosen for survival response. 
#' Probability of survival greater than \code{plot.time} is plotted as a function of cross-validated 
#' score and Treatment. See 'Details'.
#' @param perm.test Logical. If \code{perm.test=TRUE}, a permutation based test of significance
#' is conducted for statistics computed from the cross-validated scores. See 'Value' and the package 
#' vignette and for more details on the permutation tests.
#' @param nperm The number of permutations for the permutation test. Ignored if \code{perm.test=FALSE}
#' 
#' @return The return object is of class \code{eval.cv} and is a list whose components depend on the family ('cox' or 'binomial') and the 
#' chosen evaluation method ('continuous' or 'discrete')
#' @return \item{LR.Benefit}{For \code{family='cox'} and \code{method='discrete'}. The log-rank statistic for the 
#' survival difference between E and C for the 'benefit' from E group.}
#' @return \item{LR.NoBenefit}{For \code{family='cox'} and \code{method='discrete'}. The log-rank statistic for the 
#' survival difference between E and C for the 'do not benefit' from E group.}
#' 
#' @return \item{RR.T.Benefit}{For \code{family='binomial'} and \code{method='discrete'}. 
#' The response rate for subjects getting E in the 'benefit' from E group.}
#' @return \item{RR.C.Benefit}{For \code{family='binomial'} and \code{method='discrete'}. 
#' The response rate for subjects getting C in the 'benefit' from E group}
#' @return \item{RR.T.NoBenefit}{For \code{family='binomial'} and \code{method='discrete'}. 
#' The response rate for subjects getting E in the 'do not benefit' from E group}
#' @return \item{RR.C.NoBenefit}{For \code{family='binomial'} and \code{method='discrete'}. 
#' The response rate for subjects getting C in the 'do not benefit' from E group}.
#' 
#' @return \item{pval.Benefit}{If \code{perm.test=TRUE}, p-value from permutation test. 
#' For \code{family='cox'} and \code{method='discrete'}, permutation based p-value for LR.Benefit. 
#' For \code{family='binomial'} and \code{method='discrete'}, permutation based p-value for 
#' difference in response rates for E and C for the subset predicted 'benefit' from E.}
#' @return \item{pval.NoBenefit}{If \code{perm.test=TRUE}, p-value from permutation test. 
#' For \code{family='cox'} and \code{method='discrete'}, permutation based p-value for LR.NoBenefit. 
#' For \code{family='binomial'} and \code{method='discrete'}, permutation based p-value for 
#' difference in response rates for E and C for the subset predicted 'no benefit' from E.} 
#' 
#' @return \item{reg}{For \code{method='continuous'}, the regression model with treatment, predictive score and
#' treatment x predictive score interaction}
#' @return \item{pval.twosided}{For \code{method='continuous'}. Two-sided (non-directional) permutation based p-value for the treatment x predictive score
#' interaction coefficient}
#' @return \item{pval.onesided}{For \code{method='continuous'}. One-sided (directional, greater) permutation based p-value for the treatment x predictive score
#' interaction coefficient}
#' 
#' @return \item{call}{The function call}
#' 
#' @return Additional plots for both \code{method='discrete'} as well as \code{method='continuous'}. 
#' \code{print} method is available for a nice display of objects of class \code{eval.cv}. See package vignette.
#' 
#' 
#' @keywords pact, pact.cv
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' ### Survival response
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' Xf <- prostateCancer[,7:8]
#' Xv <- prostateCancer[,c(5:6,9)]
#' Treatment <- as.factor(prostateCancer[,2])
#' p <- pact.fit(Y=Y, Xf=Xf, Xv=Xv, Treatment=Treatment, family="cox", varSelect="univar")
#' cv <- pact.cv(p, nfold=5)
#' \dontrun{eval.pact.cv(cv, method="discrete", g=log(0.80), perm.test=TRUE, nperm=500)}  ## At least 20% predicted reduction in HR classified as 'sensitive'
#' eval.pact.cv(cv, method="continuous", plot.score=TRUE, perm.test=FALSE)
#' 
#' ### Binary response
#' data(EORTC10994)
#' Y <- as.factor(EORTC10994[,4])
#' Xv <- EORTC10994[,c(2,5:7)]
#' Treatment <- as.factor(EORTC10994[,3])
#' p <- pact.fit(Y=Y,Xv=Xv,Treatment=Treatment,family="binomial", varSelect="univar")
#' cv <- pact.cv(p, nfold=5)
#' \dontrun{eval.pact.cv(cv, method="discrete", g=log(1), perm.test=TRUE, nperm=500)}
#' 
#' 
### Evaluation functions for pact model

eval.pact.cv <- function(out.cv, method=c("discrete","continuous"), g=log(1), 
                         plot.score=TRUE, plot.time=NULL, perm.test=FALSE, nperm=100) {
  if(!inherits(out.cv,"pact.cv"))
    stop("First argument must be an object of class 'pact.cv'. Run pact.cv first.")
  method <- match.arg(method)
  family <- out.cv$family
  this.call <- match.call()
  if (family == "cox") {
    if (method == "discrete") {
      eval.cv <- .eval.pact.cv.survival.discrete(out.cv, g, perm.test, nperm)
    }
    else {
      if (method == "continuous") 
        eval.cv <- .eval.pact.cv.survival.continuous(out.cv, plot.score, plot.time, perm.test, nperm)
      else 
        stop("'method' can only be 'discrete' or 'continuous'")
    }
  }
  else {
    if (family == "binomial") {
      if (method == "discrete") {
        eval.cv <- .eval.pact.cv.binary.discrete(out.cv, g, perm.test, nperm)
      }
      else {
        if (method == "continuous") 
          eval.cv <- .eval.pact.cv.binary.continuous(out.cv, perm.test, nperm)
        else 
          stop("'method' can only be 'discrete' or 'continuous'")
      }
    }
    else
      stop("Error: family can only be 'cox' or 'binomial'")
  }
  eval.cv$call <- this.call
  eval.cv$method <- method
  eval.cv$family <- family
  eval.cv$perm.test <- perm.test
  eval.cv$nperm <- nperm
  class(eval.cv) <- "eval.cv"   ### Added sep 2014, needed for print method
  eval.cv
}


### Internal pact model evaluation function. For survival response. Discretized scores.
.eval.pact.cv.survival.discrete <- function(out.cv, g, perm.test, nperm) {
  predscore <- out.cv$PredScore   ### predscore represents change in log HR (E to C)
  strata <- ifelse(predscore < g, 1, 0) ### 1=(predicted) benefit from T, 0=(predicted) no benefit from E
  Y <- out.cv$Y
  SurvObj <- Surv(Y[, "time"],Y[, "status"])
  T <- out.cv$Treatment
  
  LR.Benefit <- survdiff(SurvObj ~ T, subset=(strata == 1))$chisq  ### LR statistic between control and new treatment groups for cases predicted to benefit from E
  LR.NoBenefit <- survdiff(SurvObj ~ T,subset=(strata == 0))$chisq  ### LR statistic between control and new treatment groups for cases predicted to not benefit from E
  
  par(mfrow=c(1,2), font.lab=2, cex.lab=1, cex.main=0.9)
  sf <- survfit(SurvObj ~ T, subset=(strata == 1)) ### KM plot by treatment for cases predicted to benefit from new
  plot(sf, col=c("red","blue"), lwd=2, main="Subset predicted benefit from E",
       xlab="Time", ylab="Proportion alive")
  legend("bottomleft", c("Control (C)", "New Treatment (E)"),col=c("red","blue"),lty=1,lwd=2,bty="n",cex=0.9)

  sf <- survfit(SurvObj ~ T, subset=(strata == 0)) ### KM plot by treatment for cases predicted to not benefit from new
  plot(sf, col=c("red","blue"), lwd=2, main="Subset predicted no benefit from E",
       xlab="Time", ylab="Proportion alive")
  legend("bottomleft", c("Control (C)", "New Treatment (E)"),col=c("red","blue"),lty=1,lwd=2,bty="n",cex=0.9)
  
  if (perm.test) {  ### if permutation testing of LR statistic is desired
    Xf <- out.cv$Xf
    Xv <- out.cv$Xv
    nCovarf <- out.cv$nCovarf
    nCovarv <- out.cv$nCovarv
    nfold <- out.cv$nfold
    varSelect <- out.cv$varSelect
    nsig <- out.cv$nsig
    cvfolds.varSelect <- out.cv$cvfolds.varSelect
    which.lambda <- out.cv$which.lambda
    penalty.scaling <- out.cv$penalty.scaling
    
    dimxv <- dim(Xv)
    nobs <- dimxv[1]
    
    permute.LR.Benefit <- vector("numeric", nperm)
    permute.LR.NoBenefit <- vector("numeric", nperm)
    fail.count1 <- 0
    fail.count2 <- 0
    cat("Starting Permutations..May take a few minutes to complete.. \n\n")
    for(i in (1:nperm)) {
      #### Generate a permutation of the treatment labels
      permute.T <- T[sample(nobs,nobs)]
      #### Cross-validated predictions for the permuted data
      permute.PredTmnt.CV <- .pact.cv.survival(Y, Xf, Xv, permute.T, nCovarf, nCovarv, nfold, 
                                               varSelect, nsig, cvfolds.varSelect, which.lambda, 
                                               penalty.scaling)
      permute.strata <- ifelse(permute.PredTmnt.CV$PredScore < g, 1, 0)  
      tryCatch(permute.LR.Benefit[i] <- survdiff(SurvObj ~ permute.T, subset=(permute.strata == 1))$chisq, 
               error=function(err) fail.count1 <<- fail.count1+1)
      tryCatch(permute.LR.NoBenefit[i] <- survdiff(SurvObj ~ permute.T, subset=(permute.strata == 0))$chisq, 
               error=function(err) fail.count2 <<- fail.count2+1)
    }
    cat("End Permutations \n\n")
    ### Very small or very large values of 'g' would result in too few subjects in a strata and in turn to a failed survfit model
    if (fail.count1*100/nperm > 5)  
      warning("Permutation results may not be valid for the subset predicted to benefit: g may be too small")
    pval.LR.Benefit <- (1+sum((permute.LR.Benefit > LR.Benefit), na.rm=TRUE))/(1+nperm-fail.count1)
    if (fail.count2*100/nperm > 5)
      warning("Permutation results may not be valid for the subset predicted to not benefit: g may be too large")
    pval.LR.NoBenefit <- (1+sum((permute.LR.NoBenefit > LR.NoBenefit), na.rm=TRUE))/(1+nperm-fail.count2)
    eval.cv.res <- list(LR.Benefit=LR.Benefit,LR.NoBenefit=LR.NoBenefit,
                        pval.Benefit=pval.LR.Benefit,pval.NoBenefit=pval.LR.NoBenefit)
  } else ### No permulations, no p-values
    eval.cv.res <- list(LR.Benefit=LR.Benefit,LR.NoBenefit=LR.NoBenefit)
  eval.cv.res  
}

### Internal cross-validation evaluation function. For binary response. Discretized scores.

.eval.pact.cv.binary.discrete <- function(out.cv, g, perm.test, nperm) {
  predscore <- out.cv$PredScore   ### predscore represents change in log odds (E to C)
  strata <- ifelse(predscore > g, 1, 0) ### 1=(predicted) benefit from E, 0=(predicted) no benefit from E

  Response <- out.cv$Y
  T <- out.cv$Treatment
  Sens.subset <- strata == 1
  NotSens.subset <- strata == 0
  
  #### Difference in response rate between control and new treatment groups for the 
  #### sensitive subset (subset predicted to benefit from E)
  
  Response.T.Benefit <- sum(Sens.subset & Response == 1 & T == 1, na.rm=TRUE)
  n.T.Benefit <- sum(Sens.subset & T == 1, na.rm=TRUE)

  Response.C.Benefit <- sum(Sens.subset & Response == 1 & T == 0, na.rm=TRUE)
  n.C.Benefit <- sum(Sens.subset & T == 0, na.rm=TRUE)

  test.Benefit <- prop.test(x=c(Response.T.Benefit, Response.C.Benefit), n=c(n.T.Benefit,n.C.Benefit))
  Benefit.teststat <- test.Benefit$statistic
  
  #### Difference in response rate between control and new treatment groups for the 
  #### not sensitive subset (subset predicted to not benefit from E)
  
  Response.T.NoBenefit <- sum(NotSens.subset & Response == 1 & T == 1, na.rm=TRUE)
  n.T.NoBenefit <- sum(NotSens.subset & T == 1, na.rm=TRUE)
  Response.C.NoBenefit <- sum(NotSens.subset & Response == 1 & T == 0, na.rm=TRUE)
  n.C.NoBenefit <- sum(NotSens.subset & T == 0, na.rm=TRUE)
  
  test.NoBenefit <- prop.test(x=c(Response.T.NoBenefit, Response.C.NoBenefit), n=c(n.T.NoBenefit,n.C.NoBenefit))
  NoBenefit.teststat <- test.NoBenefit$statistic 
  
  if (perm.test) {  ### if permutation testing of test statistics are desired
    Xf <- out.cv$Xf
    Xv <- out.cv$Xv
    nCovarf <- out.cv$nCovarf
    nCovarv <- out.cv$nCovarv
    nfold <- out.cv$nfold
    varSelect <- out.cv$varSelect
    nsig <- out.cv$nsig
    cvfolds.varSelect <- out.cv$cvfolds.varSelect
    which.lambda <- out.cv$which.lambda
    penalty.scaling <- out.cv$penalty.scaling
    
    dimxv <- dim(Xv)
    nobs <- dimxv[1]
    
    permute.teststat.Benefit <- vector("numeric", nperm)
    permute.teststat.NoBenefit <- vector("numeric", nperm)
    fail.count1 <- 0
    fail.count2 <- 0
    cat("Start Permutations..May take a few minutes to complete.. \n\n")
    for(i in (1:nperm)) {
      #### Generate a permutation of the treatment labels
      permute.T <- T[sample(nobs,nobs)]
      #### Cross-validated predictions for the permuted data
      permute.PredTmnt.CV <- .pact.cv.binary(Response, Xf, Xv, permute.T, nCovarf, nCovarv, nfold, 
                                             varSelect, nsig, cvfolds.varSelect, which.lambda,
                                             penalty.scaling)
      permute.strata <- ifelse(permute.PredTmnt.CV$PredScore > g, 1, 0) 
      
      permute.sens.subset <- permute.strata == 1
      permute.NotSens.subset <- permute.strata == 0
      
      Response.T.Benefit <- sum(permute.sens.subset & Response == 1 & permute.T == 1,na.rm=TRUE)
      n.T.Benefit <- sum(permute.sens.subset & permute.T == 1,na.rm=TRUE)
      Response.C.Benefit <- sum(permute.sens.subset & Response == 1 & permute.T == 0,na.rm=TRUE)
      n.C.Benefit <- sum(permute.sens.subset & permute.T == 0,na.rm=TRUE)
      
      ### tryCatch Added sep. 2014, similar to survival case
      tryCatch(permute.teststat.Benefit[i] <- prop.test(x=c(Response.T.Benefit, Response.C.Benefit), n=c(n.T.Benefit,n.C.Benefit))$statistic,
               error=function(err) fail.count1 <<- fail.count1+1) ## the value of Pearson's chi-squared test statistic.
      
      Response.T.NoBenefit <- sum(permute.NotSens.subset & Response == 1 & permute.T == 1,na.rm=TRUE)
      n.T.NoBenefit <- sum(permute.NotSens.subset & permute.T == 1,na.rm=TRUE)
      Response.C.NoBenefit <- sum(permute.NotSens.subset & Response == 1 & permute.T == 0,na.rm=TRUE)
      n.C.NoBenefit <- sum(permute.NotSens.subset & permute.T == 0,na.rm=TRUE)
      
      tryCatch(permute.teststat.NoBenefit[i] <- prop.test(x=c(Response.T.NoBenefit, Response.C.NoBenefit), n=c(n.T.NoBenefit, n.C.NoBenefit))$statistic,
               error=function(err) fail.count2 <<- fail.count2+1) ## the value of Pearson's chi-squared test statistic.
    }
    cat("End Permutations \n\n")
    ### Very small or very large values of 'g' would result in too few subjects in a strata and in turn to a failed chisq test
    if (fail.count1*100/nperm > 5)  
      warning("Permutation results may not be valid for the subset predicted to benefit: g may be too large")
    pval.Benefit <- (1+sum(permute.teststat.Benefit > Benefit.teststat))/(1+nperm)
    if (fail.count2*100/nperm > 5)
      warning("Permutation results may not be valid for the subset predicted to not benefit: g may be too small")
    pval.NoBenefit <- (1+sum(permute.teststat.NoBenefit > NoBenefit.teststat))/(1+nperm)
    eval.cv.res <- list(RR.E.Benefit=test.Benefit$estimate[1], RR.C.Benefit=test.Benefit$estimate[2],
                        RR.E.NoBenefit=test.NoBenefit$estimate[1], RR.C.NoBenefit=test.NoBenefit$estimate[2],
                        pval.Benefit=pval.Benefit, pval.NoBenefit=pval.NoBenefit)
  } else ### No permutations, no p-values
    eval.cv.res <- list(RR.E.Benefit=test.Benefit$estimate[1], RR.C.Benefit=test.Benefit$estimate[2],
                        RR.E.NoBenefit=test.NoBenefit$estimate[1], RR.C.NoBenefit=test.NoBenefit$estimate[2])
  eval.cv.res  
}


.eval.pact.cv.survival.continuous <- function(out.cv, plot.score, plot.time, perm.test, nperm) {
  predscore <- out.cv$PredScore   ### predscore represents change in log HR (E to C)
  Y <- out.cv$Y
  SurvObj <- Surv(Y[, "time"],Y[, "status"])
  T <- out.cv$Treatment
  
  CoxReg <- coxph(SurvObj ~ T + predscore + T*predscore, na.action="na.exclude")

  if ((plot.score)) { ## plot at 20%, 40%, 60% and 80% quantiles of scores
    score <- quantile(predscore, probs=c(0.2,0.4,0.6,0.8), na.rm=TRUE)
    lenq <- 4
### one plot per quantile as function of time
    newdata <- data.frame(T=as.factor(c(rep(c(1,0),lenq))), predscore=rep(score,each=2))

    mat <- matrix(c(1:4,5,5),nrow=3,ncol=2,byrow=TRUE)
    layout(mat=mat, heights = c(rep.int(0.5,2),0.1))
    par(mar=c(4, 4, 4, 2)+0.1, font.lab=2, cex.lab=1, cex.main=1, cex.axis=0.9)
    sf <- survfit(CoxReg, newdata=newdata) 
    for (i in seq(1, (2L*lenq), 2)) {
      plot(sf[c(i:(i+1))], lwd=2, main=paste("Score = ",round(newdata[i,2],2),sep=""), 
         xlab="Time", ylab="Proportion alive",col=c("blue","red"))
    }
    par(mai=c(0,0,0,0))   ### dummy plot for legend
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    legend("center",legend=c("Control (C)", "New Treatment (E)"),col=c("red", "blue"),
         lty=1,lwd=2,xpd=NA,bty="n",cex=1.1)
  }
  
  max.time <- max(Y[, "time"], na.rm=TRUE)
  if (!is.null(plot.time)){
    if (plot.time > max.time) 
      stop("Requested 'time' for plot exceeds maximum follow-up time in dataset")
    else {
      #dev.new()
      par(mfrow=c(1,1), mar=c(4, 4, 4, 5)+0.1, font.lab=2, cex.lab=0.8, cex.main=0.9, cex.axis=0.8)
      q <- seq(quantile(predscore,0.005,na.rm=TRUE), quantile(predscore,0.995,na.rm=TRUE), length.out=200)   ### generate a grid of score values
      lenq <- length(q)
      newdata <- data.frame(T=as.factor(c(rep(c(1,0),lenq))), predscore=rep(q,each=2))
      sf <- summary(survfit(CoxReg, newdata=newdata))
      prob.E <- unlist(lapply(1:lenq, function(x) {
        sfn.E <- stepfun(sf$time, c(1, sf$surv[,(2*x-1)]), right=TRUE, f=1)
        val.E <- sfn.E(plot.time)
        val.E
        }))
      prob.C <- unlist(lapply(1:lenq, function(x) {
        sfn.C <- stepfun(sf$time, c(1, sf$surv[,(2*x)]), right=TRUE, f=1)
        val.C <- sfn.C(plot.time)
        val.C
        }))
      plot(prob.E ~ q, type="l", lwd=2, col="blue", xlab="Score", 
           ylab=paste("P[Survival > ",plot.time,"]",sep=""),
           main=paste("Probability of surviving beyond landmark time"),
           ylim=c(0,1))
      lines(prob.C ~ q, lwd=2, col="red")
      legend("bottomright",legend=c("Control (C)", "New Treatment (E)"),col=c("red","blue"),
             lty=1,lwd=2,xpd=NA,bty="n",cex=1)
    }
  }

  if (perm.test) {  ### if permutation testing is desired
    Xf <- out.cv$Xf
    Xv <- out.cv$Xv
    nCovarf <- out.cv$nCovarf
    nCovarv <- out.cv$nCovarv
    nfold <- out.cv$nfold
    varSelect <- out.cv$varSelect
    nsig <- out.cv$nsig
    cvfolds.varSelect <- out.cv$cvfolds.varSelect
    which.lambda <- out.cv$which.lambda
    penalty.scaling <- out.cv$penalty.scaling
    
    dimxv <- dim(Xv)
    nobs <-  dimxv[1]
    
    intercoef <- CoxReg$coef[3] ### The true interaction coeff
    permute.intercoef <- vector("numeric", nperm)
    cat("Start Permutations..May take a few minutes to complete.. \n\n")
    for(i in (1:nperm)) {
      #### Generate a permutation of the treatment labels
      permute.T <- T[sample(nobs,nobs)]
      permute.PredTmnt.CV <- .pact.cv.survival(Y, Xf, Xv, permute.T, nCovarf, nCovarv, nfold, 
                                               varSelect, nsig, cvfolds.varSelect, which.lambda,
                                               penalty.scaling)
      permute.predscore <- permute.PredTmnt.CV$PredScore
      #### predscore model for permuted data
      permute.CoxReg <- coxph(SurvObj ~ permute.T + permute.predscore + permute.T*permute.predscore, 
                              na.action="na.exclude")
      permute.intercoef[i] <- permute.CoxReg$coef[3]
    }
    pval.twosided <- (1+sum(abs(permute.intercoef) > abs(intercoef), na.rm=TRUE))/(1+nperm)    
    pval.onesided <- (1+sum(permute.intercoef > intercoef, na.rm=TRUE))/(1+nperm)
    cat("End Permutations\n\n")
  }
  if (perm.test)  
    eval.cv.res <- list(reg=summary(CoxReg),pval.twosided=pval.twosided,pval.onesided=pval.onesided)
  else ## No permutations, no p-values
    eval.cv.res <- list(reg=summary(CoxReg))
  eval.cv.res  
}

### Internal cross-validation evaluation function. For binary response. Discretized scores.

.eval.pact.cv.binary.continuous <- function(out.cv, perm.test, nperm) {
  predscore <- out.cv$PredScore   ### predscore represents change in log odds (E to C)
  Response <- out.cv$Y
  T <- out.cv$Treatment
  
  LogReg <- glm(Response ~ T + predscore + T*predscore, 
                family = binomial(link = "logit"), na.action="na.exclude")
  
  ### Probabiliy plot
  par(mfrow=c(1,1), mar=c(4, 4, 4, 5)+0.1, font.lab=2, cex.lab=0.8, cex.main=0.9, cex.axis=0.8)
  ### generate a grid of score values
  q <- seq(quantile(predscore,0.005,na.rm=TRUE), quantile(predscore,0.995,na.rm=TRUE), length.out=200)   
  lenq <- length(q)
  newdata.E <- data.frame(T=as.factor(rep(1,lenq)), predscore=q)
  newdata.C <- data.frame(T=as.factor(rep(0,lenq)), predscore=q)
  ### Need to get prob.E and prob.C (prob of response with E and C respectively) for newdata
  prob.E <- predict(LogReg,newdata.E,type="response")
  prob.C <- predict(LogReg,newdata.C,type="response")
  plot(prob.E ~ q, type="l", lwd=2, col="blue", xlab="Score", ylab=paste("Probability of response"),
       main=paste("Probability of Response"), ylim=c(0,1))
  lines(prob.C ~ q, lwd=2, col="red")
  legend("bottomright",legend=c("Control (C)","New Treatment (E)"),col=c("red","blue"),
         lty=1,lwd=2,xpd=NA,bty="n",cex=0.9)
  
  if (perm.test) {  ### if permutation testing is desired
    Xf <- out.cv$Xf
    Xv <- out.cv$Xv
    nCovarf <- out.cv$nCovarf
    nCovarv <- out.cv$nCovarv
    nfold <- out.cv$nfold
    varSelect <- out.cv$varSelect
    nsig <- out.cv$nsig
    cvfolds.varSelect <- out.cv$cvfolds.varSelect
    which.lambda <- out.cv$which.lambda
    penalty.scaling <- out.cv$penalty.scaling
    
    dimxv <- dim(Xv)
    nobs <- dimxv[1]
    
    intercoef <- LogReg$coef[4] ### The true interaction coeff
    permute.intercoef <- vector("numeric", nperm)
    cat("Start Permutations..May take a few minutes to complete.. \n\n")
    for(i in (1:nperm)) {
      #### Generate a permutation of the treatment labels
      permute.T <- T[sample(nobs,nobs)]
      permute.PredTmnt.CV <- .pact.cv.binary(Response, Xf, Xv, permute.T, nCovarf, nCovarv, nfold, 
                                             varSelect, nsig, cvfolds.varSelect, which.lambda,
                                             penalty.scaling)
      permute.predscore <- permute.PredTmnt.CV$PredScore
      
      #### predscore model for permuted data
      permute.LogReg <- glm(Response ~ permute.T + permute.predscore + permute.T*permute.predscore, 
                            family = binomial(link = "logit"), na.action="na.exclude")
      
      permute.intercoef[i] <- permute.LogReg$coef[4]
    }
    pval.twosided <- (1+sum(abs(permute.intercoef) > abs(intercoef), na.rm=TRUE))/(1+nperm)    
    pval.onesided <- (1+sum(permute.intercoef > intercoef, na.rm=TRUE))/(1+nperm)
    cat("End Permutations \n\n")
  }
  if (perm.test)  
    eval.cv.res <- list(reg=summary(LogReg),pval.twosided=pval.twosided,pval.onesided=pval.onesided)
  else ## No permutations, no p-values
    eval.cv.res <- list(reg=summary(LogReg))
  eval.cv.res  
}


#' @title Print an object of class 'eval.cv' 
#'
#' @description
#' print method for objects of class 'eval.cv'
#' 
#' @details
#' The call that produced the object is printed, followed by the evaluation statistics. 
#'  The p-values are printed if permutation testing was asked for.
#' 
#' @method print eval.cv
#' 
#' @param x The object returned from \code{'eval.pact.cv'}
#' @param digits significant digits in the print 
#' @param ... Additional print arguments
#' 
#' @return The statistics comparing treatments E and C is printed. The printed statistics differs
#' according to whether \code{method} was \code{'discrete'} or \code{'continuous'}. p-values are 
#' printed if \code{perm.test=TRUE}
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' Xf <- prostateCancer[,7:8]
#' Xv <- prostateCancer[,c(5:6,9)]
#' Treatment <- as.factor(prostateCancer[,2])

#' p <- pact.fit(Y=Y, Xf=Xf, Xv=Xv, Treatment=Treatment, family="cox", varSelect="lasso")
#' cv <- pact.cv(p, nfold=5)
#' eval <- eval.pact.cv(cv, method="discrete", g=log(0.80))
#' eval

### Added Sep 2014
print.eval.cv <- function(x, digits = max(3, getOption("digits") - 3), ...) {  ### x = object of class 'eval.cv'
  family <- x$family
  method <- x$method
  perm.test <- x$perm.test
  nperm <- x$nperm
  
  cat("\nCall: ", deparse(x$call), "\n\n")
  cat("\nfamily: ", x$family, "\n\n")
  
  if (family == "cox") {
    if (method == "discrete") {
      cat(paste("Log-rank statistic (LR) comparing E and C in group predicted to benefit from E: ",
                round(x$LR.Benefit, digits=digits),"\n", sep=""))
      cat(paste("Log-rank statistic (LR) comparing E and C in group predicted to not benefit from E: ",
                round(x$LR.NoBenefit, digits=digits),"\n\n", sep=""))
      if(perm.test) {
        cat(paste("p-value for the LR in group predicted to benefit from E","\n","based on ",nperm,
                  " permutations: ",round(x$pval.Benefit, digits=digits),"\n\n", sep=""))
        cat(paste("p-value for the LR in group predicted to not benefit from E","\n","based on ",nperm,
                  " permutations: ",round(x$pval.NoBenefit, digits=digits),"\n\n", sep=""))
       }
    } else {  ## method="continuous"
      cat(paste("Coefficients from the regression model with Treatment, cross-validated score","\n",
          "and Treatment*score interaction","\n", sep=""))
      s <- x$reg
      print(s$coef[,1], digits=digits)
      cat("\n")
      if(perm.test) {
        cat(paste("Two-sided p-value for the Treatment*score interaction coefficient","\n",
                  "based on ",nperm," permutations: ",round(x$pval.twosided, digits=digits),"\n\n", sep=""))
        cat(paste("One-sided p-value for the Treatment*score interaction coefficient","\n",
                  "based on ",nperm," permutations: ",round(x$pval.onesided, digits=digits),"\n\n", sep=""))
      }
    }
  } else {  ## family="binomial"
    if (method == "discrete") {
      cat(paste("Response rate (RR) with E in group predicted to benefit from E: ",
                round(x$RR.E.Benefit, digits=digits),"\n", sep=""))
      cat(paste("Response rate (RR) with C in group predicted to benefit from E: ",
                round(x$RR.C.Benefit, digits=digits),"\n", sep=""))
      cat(paste("Response rate (RR) with E in group predicted no benefit from E: ",
                round(x$RR.E.NoBenefit, digits=digits),"\n", sep=""))
      cat(paste("Response rate (RR) with C in group predicted no benefit from E: ",
                round(x$RR.C.NoBenefit, digits=digits),"\n\n", sep=""))
      if(perm.test) {
        cat(paste("p-value for the difference in RR (E vs C) in group predicted to benefit from E",
                  "\n","based on ",nperm," permutations: ",round(x$pval.Benefit, digits=digits),"\n\n", sep=""))
        cat(paste("p-value for the difference in RR (E vs C) in group predicted to not benefit from E",
                  "\n","based on ",nperm," permutations: ",round(x$pval.NoBenefit, digits=digits),"\n\n", sep=""))
      }
    } else {  ## method="continuous"
      cat(paste("Coefficients from the regression model with Treatment, cross-validated score","\n",
                "and Treatment*score interaction","\n", sep=""))
      s <- x$reg
      print(s$coef[,1], digits=digits)
      cat("\n")
      if(perm.test) {
        cat(paste("Two-sided p-value for the Treatment*score interaction coefficient","\n",
            "based on ",nperm," permutations: ",round(x$pval.twosided, digits=digits),"\n\n", sep=""))
        cat(paste("One-sided p-value for the Treatment*score interaction coefficient","\n",
                  "based on ",nperm," permutations: ",round(x$pval.onesided, digits=digits),"\n\n", sep=""))
      }
    }
  }
}

### Function for overall analysis

#' @title Overall statistics and inference 
#'
#' @description
#' Produces some statistics for the overall (non-predictive) comparison of the E and C 
#' for the same dataset for which the predictive model \code{pact.fit} was developed. 
#' 
#' @details
#' Statistics for the overall comparison of the E and C is produced for the the 
#' data from a randomized clinical trial. The input is an object of class \code{pact}. 
#' 
#' @param p An object of class \code{pact}
#' 
#' @return An list with the following components. As a side effect, these are also printed
#' on screen
#' 
#' @return \item{family}{The response variable type used}
#' @return \item{nobs}{The sample size}
#' @return \item{n.E}{Number of subjects getting the new treatment}
#' @return \item{n.C}{Number of subjects getting the control treatment}
#' @return \item{LR}{The log-rank statistic for the overal difference in survival between
#' E and C groups (for family="cox")}
#' @return \item{LR.pval}{The p-value for LR based on the log-rank test (for family="cox")}
#' @return \item{RR.E}{The response rate for group treated with E (new treatment) (for family="binomial")}
#' @return \item{RR.C}{The response rate for group treated with C (Control) (for family="binomial")}
#' @return \item{RRdiff.pval}{The chi-square test based pvalue for the difference in response rates 
#' (for family="binomial")}
#' 
#' @keywords pact
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' Xf <- prostateCancer[,7:8]
#' Xv <- prostateCancer[,c(5:6,9)]
#' Treatment <- as.factor(prostateCancer[,2])
#' p <- pact.fit(Y=Y, Xf=Xf, Xv=Xv, Treatment=Treatment, family="cox", varSelect="none")
#' overall.analysis(p)
#' 
#' ### Binary response
#' data(EORTC10994)
#' Y <- as.factor(EORTC10994[,4])
#' Xv <- EORTC10994[,c(2,5:7)]
#' Treatment <- as.factor(EORTC10994[,3])
#' p <- pact.fit(Y=Y,Xv=Xv,Treatment=Treatment,family="binomial",varSelect="none")
#' overall.analysis(p)


overall.analysis <- function(p) { ### p=an object of class 'pact'  
  if(!inherits(p,"pact"))
  stop("Argument must be an object of class 'pact'.")
  
  Y <- p$Y
  T <- p$Treatment
  family <- p$family
  nobs <- dim(as.matrix(Y))[1]
  
  n.T <- sum(T == 1,na.rm=TRUE)
  n.C <- sum(T == 0,na.rm=TRUE)
  
  if (family == "binomial") {
    Response.T <- sum(Y == 1 & T == 1,na.rm=TRUE)
    Response.C <- sum(Y == 1 & T == 0,na.rm=TRUE)
    Overall.test <- prop.test(x=c(Response.T, Response.C), n=c(n.T,n.C))
    RR.T <- Overall.test$estimate[1]
    RR.C <- Overall.test$estimate[2]
    pval <- Overall.test$p.value
    output <- paste("Description of the problem:\n",
                    "--------------------------\n\n",
                    "family: ", family,
                    "\n\nNumber of subjects: ", nobs,
                    "\n\nNumber of subjects assigned new treatment: ", n.T,
                    "\n\nNumber of subjects assigned standard treatment: ", n.C,
                    "\n\nOverall Analysis \n",
                    "--------------------\n\n",
                    "Response rate in the group administered new treatment: ",
                    round(RR.T,2),
                    "\n\nResponse rate in the group administered standard treatment: ",
                    round(RR.C,2),
                    "\n\n(Chi-square test based) P-value for the test of significance of the difference in response rates: ",
                    round(pval,4),
                    "\n\n",
                    sep="")
    cat(output)
    outlist <- list(family=family, nobs=nobs, n.E=n.T, n.C=n.C, RR.E=RR.T, RR.C=RR.C, RRdiff.pval=pval)
  } else {
    SurvObj <- Surv(Y[, "time"],Y[, "status"])
    sf <- survfit(SurvObj ~ Treatment)
    LR.full <- survdiff(SurvObj ~ Treatment)$chisq
    LR.pval <- 1 - pchisq(LR.full, 1)
    
    output <- paste("\nDescription of the problem:\n",
                    "----------------------------\n\n",
                    "family: ", family,
                    "\n\nNumber of subjects: ", nobs,
                    "\n\nNumber of subjects assigned new treatment: ", n.T,
                    "\n\nNumber of subjects assigned standard treatment: ", n.C,
                    "\n\nOverall Analysis:\n",
                    "--------------------\n\n",
                    "Log-rank statistic: ", round(LR.full,3),
                    "\n\np-value (Log-rank test): ", round(LR.pval,2),
                    "\n\nOverall Kaplan-Meier plot in the plot window \n\n",sep="")
    cat(output)
    
    par(font.lab=2, cex.lab=1.2)
    plot(sf, col=c("red","blue"), lwd=2, main="Overall Analysis",
         xlab="Time", ylab="Proportion alive")
    legend("bottomleft", c("Standard Treatment", "New Treatment"), col=c("red","blue"), lty=1, lwd=2, bty="n", cex=1)
    outlist <- list(family=family, nobs=nobs, n.E=n.T, n.C=n.C, LR=LR.full, LR.pval=LR.pval)
  }
  ##outlist
}

#' Prostate cancer dataset
#' 
#' A dataset containing survival information for 485 subjects with prostate cancer. 
#' See 'Details' for the variables in this dataset.
#' 
#' \itemize{
#'   \item ID. Subject identifier 
#'   \item Treatment. Treatment received. '0' for control and '1' for new
#'   \item time. Survival time in months
#'   \item status. Censoring status. '0' for censored, '1' for died
#'   \item age. Age (in years) at diagnosis
#'   \item pf. Performance status (Normal Activity or Limited Activity)
#'   \item sz. Size of the primary tumor (cm2)
#'   \item sg. Index of a combination of tumor stage and histologic grade 
#'   \item ap. Serum phosphatic acid phosphatase levels
#'   }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 485 rows and 9 variables
#' @name prostateCancer
NULL

#' EORTC10994 dataset
#' 
#' A dataset containing treatment, response and covariate information for 125 subjects 
#' with breast cancer. See 'Details' for the variables in this dataset.
#' 
#' \itemize{
#'   \item ID. Subject identifier 
#'   \item Treatment. Treatment received. '0' for control and '1' for experimental
#'   \item Response. Binary response to treatment - '0' for 'non-responder' and '1' for 'responder'
#'   \item Age. Age (in years) at diagnosis
#'   \item TumorSize. size of tumor
#'   \item Node. Node positive (Yes) or node negative (No)
#'   \item ERBB2Log2. log2 transformed ERBB2 levels
#'   }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 125 rows and 7 variables
#' @name EORTC10994
NULL

#' GSE10846 dataset
#' 
#' A dataset containing the survival, treatment and gene expression data for 
#' 412 patients with diffuse large B cell lymphoma and treated with CHOP or 
#' CHOP+Rituximab
#' 
#' \itemize{
#'   \item time. Survival time in months
#'   \item status. Censoring status. '0' for censored, '1' for died
#'   \item Treatment. Treatment received. '0' for control (CHOP) and '1' 
#'   for new (CHOP+Rituximab)
#'   \item Columns 4-1003 are gene expression values (normalized using
#'    the MAS 5.0 software and log2 transformed) of the first 1000 genes  
#'    highest variance from the original dataset from GEO
#'   \item The row names are the array names
#'   }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 412 rows and 1003 columns
#' @name GSE10846
NULL