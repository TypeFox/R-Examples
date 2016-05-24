####### Here the functions related to calculating lsmeans and difflsmeans ######
####### and for calculating Satterthwaiteare are specified #####################
####### these are the modifications of lmerTest functions ######################
####### the modified versions are needed to be able do post-hoc for MAM ########
####### see explanations in paper and thesis ###################################



getX <- function(model)
{
  return(getME(model, "X"))
}

rhoInitLsmeans <- function(model, modelMAM = NULL){
  # creating rho
  rho <- new.env(parent = emptyenv()) # create an empty environment
  rho$model <- model
  if(!is.null(modelMAM)){
    rho$modelMAM <- modelMAM
    rho$sigma <- sigma(modelMAM)
    rho$thopt <- getME(modelMAM, "theta")
    ## save optima from the MAM
    rho$theopt <- c(rho$thopt, rho$sigma)
    dd <- devfun5(model,  getME(model, "is_REML"))
    h <- myhess(dd, rho$theopt)
    ## calculate variance-covariance matrix of random effects
    rho$A <- 2*solve(h)
  }
  else{
    rho$sigma <- sigma(model)
    rho$thopt <- getME(model, "theta")
    ## save optima from the MAM
    rho$theopt <- c(rho$thopt, rho$sigma)
    dd <- devfun5(model,  getME(model, "is_REML"))
    h <- myhess(dd, rho$theopt)
    ## calculate variance-covariance matrix of random effects
    rho$A <- 2*solve(h)
  }
    
  rho$X <- getX(model) 
  
  rho$REML <-  getREML(model)  

  rho$fixEffs <- fixef(model)  

  return(rho)  
}


###############################################################################
# function to calculate T test JSS
###############################################################################
calculateTtestJSS <- function(rho, Lc, nrow.res, ddf="Satterthwaite")
{
  
  resultTtest <- matrix(0, nrow = nrow.res, ncol = 4)
  colnames(resultTtest) <- c("df", "t value", "p-value", "sqrt.varcor")
  
  if(ddf == "Kenward-Roger"){
    if (!requireNamespace("pbkrtest", quitly = TRUE)) 
      stop("pbkrtest package required for Kenward-Roger's approximations")
    Va <- pbkrtest::vcovAdj(rho$model)
  }
  else{
    # based on theta parameters
    vss <- vcovJSStheta2(rho$model)
    # based on variance parameters
    #vss <- vcovJSStheta2.var(rho$model)
  }
  
  for(i in 1:nrow.res)
  {
    
    if(ddf == "Kenward-Roger"){
      L <- Lc[,i]
      .ddf <- pbkrtest::get_ddf_Lb(rho$model, L)      
      b.hat <- rho$fixEffs
      Lb.hat <- sum(L * b.hat)
      Va.Lb.hat <- t(L) %*% Va %*% L
      t.stat <- as.numeric(Lb.hat / sqrt(Va.Lb.hat))
      p.value <- 2 * pt(abs(t.stat), df = .ddf, lower.tail = FALSE)
      resultTtest[i,1] <- .ddf
      resultTtest[i,2] <- t.stat
      resultTtest[i,3] <- p.value
      resultTtest[i,4] <- as.numeric(sqrt(Va.Lb.hat))
    }
    else{
      ## based on theta parameters
      g <- mygrad(function(x)  vss(t(Lc[,i]), x), c(rho$thopt, rho$sigma))
      ## based on var cor parameters
      #g <- grad(function(x)  vss(t(Lc[,i]), x), rho$vars)
      
      #denominator df
      denom <- t(g) %*% rho$A %*% g
      ## for the theta and sigma parameters
      varcor <- vss(t(Lc[,i]), c(rho$thopt, rho$sigma))
      ## for var cor parameters
      #varcor <- vss(t(Lc[,i]), rho$vars)
      #df
      resultTtest[i,1] <- 2*(varcor)^2/denom
      #statistics
      resultTtest[i,2] <- (Lc[,i] %*%rho$fixEffs)/sqrt(varcor) 
      resultTtest[i,3] <- 2*(1 - pt(abs(resultTtest[i,2]), df = resultTtest[i,1]))
      resultTtest[i,4] <- sqrt(varcor) 
    }
    
  }
  
  return(resultTtest)
}

## from lmerTest package version 2.0-24
## with deleted update of contrasts
refitLM <- function(obj) {
  
  mm <- model.frame(obj)
  colnames(mm)[1] <- "y"
  fo <- getFormula(obj, withRand=FALSE)# formula(obj,fixed.only=TRUE)
  fo <- update(fo, y ~ .)
  lm(fo, data=mm)
}

############################################################################
#get formula for model: function from lmerTest version 2.0-24
############################################################################
getFormula <- function(model, withRand=TRUE)
{
  fmodel <- formula(model)  
  
  if(withRand)
    return(fmodel)
  
  fm <- paste(fmodel)
  fmodel.red <- paste(fm[2],fm[1], 
                      paste(fm[3], 
                            paste(unlist(lapply(names(.fixedrand(model)$randeffs),
                                                function(x) paste("(",x, ")"))), 
                                  collapse = " - "), 
                            sep = "-"))
  return(update(fmodel, fmodel.red))
}

## function from lmerTest
## the same function is in SensMixed package
## list of random and fixed terms of the model
.fixedrand <- function(model)
{  
  effs <- attr(terms(formula(model)), "term.labels")
  neffs <- length(effs)
  randeffs <- effs[grep(" | ", effs)]
  randeffs <- sapply(randeffs, function(x) substring(x, 5, nchar(x)))
  fixedeffs <- effs[!(effs %in% names(randeffs))]  
  return(list(randeffs=randeffs, fixedeffs=fixedeffs))
}