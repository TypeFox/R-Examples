### The function to fit the log multiplicative association model

### The function 'llla'
### fits the log linear by linear association model
### using pseudolikelihood method
### Syntax
###   llla(data, item.mtx, trait.mtx, useMLE)
###
### Arguments
###  'data': is a data frame or matrix with rows indicating
###   individuals and columns indicating items and the values indicating
###   the choices.
###  'item.mtx': is the adjacency matrix between items and
###   the latent traits.
###  'trait.mtx': is the adjacency matrix for latent traits
###  'useMLE': inidicates whether maximum likelihood estimation is used
###  'uncorrected': if the value is TRUE, calculate the uncorrected
###     standard errors
###
### Values:
###  'coefficients': the parameter estimates in the LLLA model
###  'se': the standard error of coefficient esimates(sandwich estimator)
###  'covb': the covariance matrix of the coefficient esimates
###  'se.uncorrected': the standard error not corrected
###  'ncat' : number of categories
###  'nitem' : number of items
###  'nexaminee' : number of examinees (may be change it to 'nperson'?)
### Details:
###   This is the main function to fit the LLLA model
###
### Author: Zhushan "Mandy" Li
### Last Updated: Sept 16, 2006
###
### Update Notes:
###  Aug 23: added maximum likelihood estimates 

llla <- function(data, item.mtx=rep(1,ncol(data)), trait.mtx=1, useMLE=FALSE, uncorrected=FALSE){
  call <- match.call();
  
  if(useMLE){ # maximum likelihood
    result <- llla.mle(data, item.mtx, trait.mtx)
  }
  else
  { # pseudolikelihood
  ## check whether responses are binary or multicategory
  ncat <- length(unique(as.matrix(data[!is.na(data)])))

  if(ncat == 2){ ## if it is binary
    result <- llla.bi(data, item.mtx, trait.mtx, uncorrected) # fit the binary model
  }
  else if(ncat > 2) { ## if it multicategory
    result <- llla.mc(data, item.mtx, trait.mtx, uncorrected) # fit the multicategory model
  } else {
    stop("Error, the data has less than 2 categories")
  }
  
  }

  result$call <- call
  result$item.mtx <- item.mtx
  result$trait.mtx <- trait.mtx
  if(useMLE){
    result$method <- "MLE"
  } else {
    result$method <- "PLE"
  }
  
  class(result) <- "llla"
  return(result)
}


#### 'llla.bi' : fitting the binary model
### Details
###   This function produce Pseudolikelihood esimates of the model for binary
###  Responses. It stacks the data first, then fits a logistic regression model
###  to the stacked data. Sandwich estimator of se's are calculated.
llla.bi <- function(data, item.mtx, trait.mtx, uncorrected=FALSE){
  ### identify the model
  ncat <- 2
  nitem <- ncol(data)
  nexaminee <- nrow(data)
  
  ### stack the data
  pldata <- plStackData(data, item.mtx, trait.mtx)
  pldata <- data.frame(pldata)

  ### extract the variable names that contain "phi"
  pldata.names <- names(pldata)
  phi.names <- pldata.names[grep("phi", pldata.names)]
  
  ### fit the logistic regression on the stacked data
  lrFormula <- as.formula( paste("resp ~ factor(item) +",
                                 paste(phi.names, collapse='+'), "-1") ) 
  plfit <- glm(lrFormula, family=binomial(logit), data=pldata)

  ### get sandwich estimator
  covb.sandwich <- vcov.sandwich.glm(plfit, strata=pldata$cid)
  se.sandwich <- sqrt(diag(covb.sandwich))

  ### get uncorrected estimator
  covb.uncorrected <- NULL;
  se.uncorrected <- NULL;
  if (uncorrected) {
  covb.uncorrected <- vcov(plfit, strata=pldata$cid)
  se.uncorrected <- sqrt(diag(covb.uncorrected))
  }
  
  ### return the results
  coefficients <- coef(plfit)

  ### rename the coefficients
  coef.item.names <- paste("item", 1:nitem, sep='')
  coef.names <- c(coef.item.names, phi.names)

  names(coefficients) <- coef.names;
  names(se.sandwich) <- coef.names;
  rownames(covb.sandwich) <- coef.names;
  colnames(covb.sandwich) <- coef.names;
  if (uncorrected){
    names(se.uncorrected) <- coef.names;
    rownames(covb.sandwich) <- coef.names;
    colnames(covb.sandwich) <- coef.names;
  }
  
  fit <- list(coefficients=coefficients, se=se.sandwich, covb=covb.sandwich,
              covb.uncorrected=covb.uncorrected, se.uncorrected=se.uncorrected,
              ncat=ncat, nitem=nitem, nexaminee=nexaminee)
  return(fit)
}


#### 'llla.mc' : fitting the multicategory model
### Details
###    This function produce Pseudolikelihood esimates of the model for multi-
###  categorical  Responses.
###    It stacks the data first, then fits a multinomial condtional logit  model
###  to the stacked data by using cox regression (function 'coxph' in package
###  'survival'.
###    Sandwich estimator of se's are given by 'coxph'.
# require('survival')

llla.mc <- function(data, item.mtx, trait.mtx, uncorrected=FALSE){
  ### identify the model
  ncat <- length(unique(as.vector(as.matrix(data))))
  nitem <- ncol(data)
  nexaminee <- nrow(data)
  
  ### stack the data
  pldata <- plStackData(data, item.mtx, trait.mtx)
  pldata <- data.frame(pldata)
  pldata$resp <- as.factor(pldata$resp)
  pldata$item <- as.factor(pldata$item)
  
  ### extract the variable names that contain "phi"
  pldata.names <- names(pldata)
  phi.names <- pldata.names[grep("phi", pldata.names)]
  nphi <- length(phi.names)
  
  ### fit the conditional logit model on the stacked data
  pc <- mclgen(pldata, 'resp')
##  require(survival)
  resp.X <- model.matrix(~pc$resp)
  resp.X<-resp.X[,attributes(resp.X)$assign==1]
  
  pc[,phi.names] <- rep(0:(ncat-1), nrow(pc)/ncat)*pc[,phi.names]

  clFormula <- as.formula(
     paste('Surv(rep(1, nrow(pc)), depvar) ~  item:resp.X + ', phi.names,
           '+strata(id)+cluster(cid)')
                          )

  plfit <- coxph(clFormula, data=pc, robust=TRUE)

  ### return the results
  covb <- vcov(plfit)
  se <- sqrt(diag(covb))
  coefficients <- coef(plfit)

  ## reorder and rename the coefficients so that they are consistent with
  ## binary case.
  p <- length(coefficients)
  reorder <- c( (nphi+1):p, 1:nphi )
  coefficients <- coefficients[reorder]
  se <- se[reorder]
  covb <- covb[reorder, reorder]

  coef.item.names <- paste("item", outer(1:nitem, 1:(ncat-1), paste, sep='.'), sep='')
  coef.names <- c(coef.item.names, phi.names)
  names(coefficients) <- coef.names
  names(se) <- coef.names
  rownames(covb) <- coef.names; colnames(covb) <- coef.names

  ### get uncorrected estimator
  covb.uncorrected <- NULL;
  se.uncorrected <- NULL;
  if (uncorrected) {
    plfit.uncorrected <- coxph(clFormula, data=pc, robust=FALSE)
    
    covb.uncorrected <- vcov(plfit.uncorrected)
    covb.uncorrected <- covb.uncorrected[reorder, reorder]
    se.uncorrected <- sqrt(diag(covb.uncorrected))

    names(se.uncorrected) <- coef.names
    rownames(covb.uncorrected) <- coef.names; colnames(covb.uncorrected) <- coef.names
  }
  
  ##
  fit <- list(coefficients=coefficients, se=se, covb=covb,
              se.uncorrected=se.uncorrected, covb.uncorrected=covb.uncorrected,
              ncat=ncat, nitem=nitem, nexaminee=nexaminee)
  return(fit)
}

#### mclgen, adapted from package 'catspec'
#### generate data set for fitting multinomial logit model with conditonal
#### logit methods by using funcion 'coxph'
mclgen <- function (datamat, catvar) 
{
    stopifnot(is.data.frame(datamat))
#    attach(datamat)
    stopifnot(is.factor(datamat[,catvar]))
    ncat <- nlevels(datamat[,catvar])
#    perschoice <- as.data.frame(rep(datamat, ncat))
    perschoice <- reshape(datamat, direction = "long", varying = lapply(names(datamat), 
        rep, ncat), timevar = "newy")
    perschoice <- perschoice[sort.list(perschoice$id), ]
    dep <- parse(text = paste("perschoice$", catvar))
    perschoice$depvar <- ifelse(as.numeric(eval(dep)) == perschoice$newy, 
        1, 0)
    perschoice[[catvar]] <- as.factor(perschoice$newy)
    perschoice[[catvar]] <- factor(eval(dep), labels = levels(datamat[,catvar]))
#    detach(datamat)
    perschoice
}


#### fitting LLLA using maximum likelihood
llla.mle <- function(data, item.mtx, trait.mtx)
{
  ##
  ncat <- length(unique(as.matrix(data[!is.na(data)])))
  nitem <- ncol(data)
  nexaminee <- nrow(data)

  ## generate the data used for MLE
  mldata <- mlData(data, item.mtx, trait.mtx)

  ## extract the variable names to be used for formula construction later
  mldata.names <- names(mldata)
  phi.idx <- grep("phi", mldata.names)
  phi.names <- mldata.names[phi.idx]
  nphi <- length(phi.idx)
  
  count.idx <- grep("count", mldata.names)
  count.name <- mldata.names[count.idx]

  item.idx <- (1:length(mldata.names))[-c(phi.idx, count.idx)] 
  item.names <- mldata.names[item.idx]
  
  ## fit the poisson regression model
  mlFormula <-  as.formula(
       paste("count ~",
             paste(paste("factor(", item.names, ")",sep=''), collapse="+"),
             "+", 
             paste(phi.names, collapse='+')
             )
                           )
  mlfit <- glm(mlFormula, family=poisson(log), data=mldata)

  ## return the results
  coefficients <- coef(mlfit)
  covb <- vcov(mlfit)
  se <- sqrt(diag(covb))

  ## reorder the coefficients, so that the order of the coefficients is
  ## the same as the order given by llla.bi and llla.mc 
  neworder <- c(1, t(matrix(2:(nitem*(ncat-1)+1), ncol=nitem)),nitem*(ncat-1)+1+(1:nphi))
  coefficients <- coefficients[neworder]
  covb <- covb[neworder, neworder]
  se <- se[neworder]
  
  ### rename the coefficients
  if(ncat==2){ # binary response
    coef.item.names <- paste("item", 1:nitem, sep='')
  } else { # multicategory repsonse
    coef.item.names <- paste("item", outer(1:nitem, 1:(ncat-1), paste, sep='.'), sep='')
  }
  
  coef.names <- c("(Intercept)",coef.item.names, phi.names)

  names(coefficients) <- coef.names;
  names(se) <- coef.names;
  rownames(covb) <- coef.names;
  colnames(covb) <- coef.names;

  ##
  fit <- list(coefficients=coefficients, se=se, covb=covb)
  return(fit)  
}
    

##### How to print the llla object
print.llla <- function(x, ...){
  cat("Log linear by linear Model\n")
  cat("Method: ",
      ifelse(x$method=="PLE", "Pseudo-likelihood Method",
             "Maximum Likelihood Method"),
      "\n\n")
  cat("Call: ",  deparse(x$call), "\n\n")
  cat("Coefficients:\n")
  print(round(cbind(Estimate=x$coefficients, Std.Error=x$se), 4))
  invisible(x)
}

summary.llla <- function(object, ...){
  return(object, ...)
}

vcov.llla <- function(object, ...){
  return(object$covb, ...)
}

coef.llla <- function(object, ...){
  return(object$coefficients)
}
