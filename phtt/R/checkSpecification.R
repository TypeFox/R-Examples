########################################################################################
## Model-Test (Kneip Sickles Song 2009 Section 3.4)  and Hausman type test of Bai (2009)
##
## obj1: Eup- or KSS-object
## obj2: Eup-object 
## level      = numeric, significance level
##
## Gives:
## Test-Statistic, p.value, crit.value and sig.level
##
########################################################################################

## rm(list=ls())

checkSpecif <- function(obj1, obj2, level = 0.05){
  if(missing(obj2)){
    test.int <- TRUE
  }else{
    test.int <- FALSE
  }
  ##===================================================================================
  ## check arguments

  ## ==================================================================================
  ## effects
  if(!test.int){
    if(class(obj1)!="Eup" & !class(obj2)!="Eup"){
      stop("\n >>obj1<< and >>obj2<< have to be Eup-objects.")
    }
    testname <- "Hausman type test as proposed by Bai (2009)"
    T  <- obj1$dat.dim[1]; if(obj2$dat.dim[1]!=T) stop("The first and the second object don't have the same data dimensions.")
    n  <- obj1$dat.dim[2]; if(obj2$dat.dim[2]!=n) stop("The first and the second object don't have the same data dimensions.")	
    P  <- obj1$dat.dim[3]; if(obj2$dat.dim[3]!=P) stop("The first and the second object don't have the same data dimensions.")
    d1 <- obj1$used.dim
    d2 <- obj2$used.dim
    
    if(d1 == 0){
      if(d2 ==0){stop("Both objects have only additive unobserved effects.")}
      add.Obj <- obj1
      int.Obj <- obj2
      additive.effects <- add.Obj$additive.effects
      if(!additive.effects %in% c("individual", "time", "twoways")) stop("'obj1' does not have classical individual, time, or twoways-effects.")
      if(d2==1 && additive.effects == "twoways") stop("The factor dimension in the model under the alternative hypothesis is wrong.\n 'factor.dim' should be 2.")
      if(d2==2 && additive.effects %in% c( "individual", "time")) stop("The factor dimension in the model under the alternative hypothesis is wrong.\n 'factor.dim' should be 1.")
      if(d2>2) stop("The factor dimension in the model under the alternative hypothesis is larger than 2.")
    }else{
      if(d2 != 0) stop("Both models have time-varying interactive effects.")
      add.Obj <- obj2
      int.Obj <- obj1
      additive.effects <- add.Obj$additive.effects
      if(!additive.effects %in% c( "individual", "time", "twoways")) stop("'obj2' does not have classical individual, time, or twoways-effects.")
      if(d1==1 && additive.effects == "twoways") stop("The factor dimension in the model under the alternative hypothesis is wrong.\n 'factor.dim' should be 2.")
      if(d1==2 && additive.effects %in% c( "individual", "time")) stop("The factor dimension in the model under the alternative hypothesis is wrong\n 'factor.dim' should be 1.")
      if(d1>2) stop("The used factor dimension in the model under the alternative hypothesis is larger than 2.")
    }
    beta.add <- add.Obj$slope.para
    C <- Eup.inference(add.Obj, error.type = 1, kernel.weights = NULL)$inv.ZZ
	
    beta.int <- int.Obj$slope.para
    infbetaint <- Eup.inference(int.Obj, error.type = 1, kernel.weights = NULL)
    D <- infbetaint$inv.ZZ
    sig2.hat <- infbetaint$sig2.hat
    DmC <- D-C
    
    
    Test.Stat <- n*T*sig2.hat^{-1}*t(beta.int - beta.add)%*%solve(DmC)%*%(beta.int - beta.add)
    if(Test.Stat < 0){stop(paste("The assumptions of the test are not fulfilled.\n The (unobserved) true number of factors is probably greater than ", max(d1,d2),".", sep=""))}
    p.value      <- 1 - pchisq(Test.Stat, df = P)
    crit.value   <- qchisq(level, df = P)
    result       <- list(Test.Stat= round(Test.Stat,2), p.value= round(p.value, 2),
                     crit.value = round(crit.value, 2), sig.level = round(level, 2))
    result$print <- additive.effects
  }else{
    if(class(obj1)=="Eup"){
      testname     <- "Test of Kneip, Sickles, and Song (2012)"
      obj          <- Eup(formula=obj1$formula, additive.effects=obj1$additive.effects, factor.dim=0)
      resObj       <- obj$unob.fact.stru + obj$residuals
      fsvd.pca.obj <- fsvd.pca(resObj, spar = 0)
      result       <- KSS.OptDim(fsvd.pca.obj, criteria = c("KSS.C"), alpha    = level)[[2]]
    }
    if(class(obj1)=="KSS"){
      testname     <- "Test of Kneip, Sickles, and Song (2012)"
      obj          <- KSS(formula=obj1$formula, additive.effects=obj1$additive.effects, factor.dim=0)
      resObj       <- obj$unob.fact.stru + obj$residuals
      fsvd.pca.obj <- fsvd.pca(resObj, spar = 0)
      result       <- KSS.OptDim(fsvd.pca.obj, criteria = c("KSS.C"), alpha    = level)[[2]]
    }
  }
  result$test.int <- test.int
  result$testname <- testname
  class(result)  <- "checkSpecif" 
  return(result)
}


## Methods ========================================================================================

print.checkSpecif <- function(x,...){
  cat("----------------------------------------------\n")
  if(!x$test.int) cat("Testing Additive vs. Interactive Effects\n")
  else  cat("Testing the Presence of Interactive Effects\n")
  cat(paste(x$testname))
  cat("\n----------------------------------------------\n")
  if(!x$test.int) cat(paste("H0: There are only additive ",x$print,"-effects\n\n",sep=""))
  else  cat(paste("H0: The factor dimension is equal to 0.\n\n",sep=""))
  outp        <- c(x$Test.Stat, signif(as.numeric(x$p.value),digits=3), x$crit.value, x$sig.level)
  names(outp) <- c("Test-Statistic", "p-value", "crit.-value", "sig.-level")
  print(outp)
}


