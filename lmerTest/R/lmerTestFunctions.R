## all functions (step, anova, rand, lsmeans, difflsmeans) 
## call first this function
## so everything starts from here
totalAnovaRandLsmeans <- function(model, ddf = "Satterthwaite", type = 3, 
                                  alpha.random = 0.1, alpha.fixed = 0.05, 
                                  reduce.fixed = TRUE, reduce.random = TRUE, 
                                  fixed.calc = TRUE, lsmeans.calc = TRUE, 
                                  difflsmeans.calc = TRUE,  isTotal = FALSE, 
                                  isAnova = FALSE, isRand = FALSE, 
                                  isLSMEANS = FALSE, 
                                  isDiffLSMEANS = FALSE, isTtest = FALSE, 
                                  test.effs = NULL, keep.effs = NULL)
{
  
  change.contr <- TRUE
  
  ## check type of hypothesis
  if(!isRand && !(type %in% c(1,2,3)))  
    stop('Parameter type is wrongly specified') 
  
  ## check keep.effs 
  if(!isTotal)
    keep.effs <- NULL
  else{    
    if(!is.null(keep.effs)){
      model.effs <- .fixedrand(model)
      keep.effs1 <- .getKeepEffs(keep.effs, model.effs) 
      if(length(unlist(keep.effs1)) == 0)
        message(paste("No ", keep.effs, "exist among effects in the model"))
      keep.effs <- keep.effs1
    }    
  }
  
  
  data <- model.frame(model) 
  
  #update contrasts for anova or step methods
  mm <- model.matrix(model)
  l.lmerTest.private.contrast<- attr(mm,"contrasts")
  contr <- l.lmerTest.private.contrast
  
  ## THE  CONTRASTS CODE IS TRANSFERRED AFTER THE REDUCTION
  ## OF THE RAND EFFECTS - THE CHANGE OF THE CONTRASTS INFLUENCED 
  ## THE LRT FOR RANDOM EFFECTS - EXAMPLE IN testContrasts.R

  #not to show the warnings  
  #options(warn=-1) 
  
  result <- NULL
  anova.table <- NULL
  
  result$response <- rownames(attr(terms(model),"factors"))[1]
  
  
   ## deleted because use ML for anova(m1, m2) for random effects
  if( isRand || isTotal || (ddf=="Kenward-Roger" && (isTotal || isAnova)) )
  {
    
      model<-
      if (getREML(model) == 1)
      {
        model
      }
      else
      {
        warning("\n model has been refitted with REML=TRUE \n")
        updateModel(model, .~., reml.lmerTest.private=TRUE, 
                    l.lmerTest.private.contrast)
      }
  }
  
  mf.final <- update.formula(formula(model),formula(model)) 
  
  # save the call of the model              
  result$call <- model@call
  
  #check if there are correlations between intercept and slope
  result$corr.intsl <- checkCorr(model)  
  
  
   
  #save results for fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
  {
    result <- saveResultsFixModel(result, model)
    result$rand.table=NULL
    return(result)
  }
  
  
  

  #analysis of the random part  
  if(isRand || isTotal)
  {
    if(isRand)
      reduce.random <- FALSE
    result.rand <- elimRandEffs(model, data, alpha.random, reduce.random, 
                                l.lmerTest.private.contrast, keep.effs$randeffs)  
   
    model <- result.rand$model
    #convert rand table to data frame
    rt <- as.data.frame(result.rand$TAB.rand)
    rt$Chi.DF <- as.integer(rt$Chi.DF)
    if(!is.null(rt$elim.num))
      rt$elim.num <- as.integer(rt$elim.num)
    
    result$rand.table <- rt
    if(isRand || !fixed.calc){
      result$model <- model
      return(result)
    }
      
  }
      
  #save results for fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
    return(saveResultsFixModel(result, model, type))


## remove update the contrasts
if(change.contr){
  ### change contrasts for F tests calculations
  #list of contrasts for factors
  if( isAnova || isTotal )
  {    
    if( length(which(unlist(contr)!="contr.SAS")) > 0 )
    {
      names.facs <- names(contr)
      l.lmerTest.private.contrast <- as.list(rep("contr.SAS",length(names.facs)))
      names(l.lmerTest.private.contrast) <- names(contr)
      model <- updateModel(model, .~., getREML(model), l.lmerTest.private.contrast) 
    }    
  }
  else
  {
    #update model to mer class
    model <- updateModel(model, .~., getREML(model), l.lmerTest.private.contrast)
  }
}
  
  
  #perform reduction of fixed effects for model with mixed effects
  stop = FALSE
  is.first.anova <- TRUE
  is.first.sign <- TRUE  
  
  
  
  while(!stop)
  {      
    
    # if there are no fixed terms
    if(nrow(anova(model, ddf="lme4"))==0)
    {
      if(is.null(anova.table))
      {
        if(isLSMEANS || isDiffLSMEANS)
        {
          lsmeans.summ <-  matrix(ncol=7,nrow=0)
          colnames(lsmeans.summ) <- c("Estimate", "Standard Error", "DF", 
                                      "t-value", "Lower CI", "Upper CI", "p-value")
          lsmeans.summ <- as.data.frame(lsmeans.summ)
          if(isLSMEANS)
            result$lsmeans.table <- lsmeans.summ
          if(isDiffLSMEANS)
            result$diffs.lsmeans.table <- lsmeans.summ
          return(result)
        }
        if(isTtest)
        {
          # save lmer outcome in rho environmental variable
          rho <- rhoInitJSS(model)  
          
          # calculate asymptotic covariance matrix A
          dd <- devfun5(model,  getME(model, "is_REML"))
          h <- myhess(dd, c(rho$thopt, sigma = rho$sigma))
          
          rho$A <- 2*solve(h)
          
                    
          tsummary <- calculateTtestJSS(rho, diag(rep(1,length(rho$fixEffs))), 
                                        length(rho$fixEffs), ddf = ddf)
          
          result$ttest <- list(df=tsummary[, "df"], tvalue=tsummary[, "t value"], 
                               tpvalue=tsummary[, "p-value"])
        }
        result$model <- model
        result$anova.table <- anova(model, ddf="lme4")
        return(result)        
        
      }          
      break
    }        
    
    
    # save lmer outcome in rho environmental variable
    rho <- rhoInitJSS(model)
  
    
    # calculate asymptotic covariance matrix A??
    if(!(ddf == "Kenward-Roger" && isAnova)){
  
    
    
      ## based on theta parameters and sigma
      # also correct
      dd <- devfun5(model,  getME(model, "is_REML"))
      h <- myhess(dd, c(rho$thopt, sigma = rho$sigma))  
      
           
      
      ch <- try(chol(h), silent=TRUE)
      if(inherits(ch, "try-error")) {
        message("Model is not identifiable...")
      }
      rho$A <- 2*chol2inv(ch)
      
      eigval <- eigen(h, symmetric=TRUE, only.values=TRUE)$values
      isposA <- TRUE
      if(min(eigval) < sqrt(.Machine$double.eps)) ## tol ~ sqrt(.Machine$double.eps)
        isposA <- FALSE
      
       
      if(!isposA)
      {
        print("Asymptotic covariance matrix A is not positive!")
      }
      
      
    }
    
    
    #calculate ttest and p-values for summary
    if(isTtest)
    {
      
      tsummary <- calculateTtestJSS(rho, diag(rep(1,length(rho$fixEffs))), 
                                  length(rho$fixEffs), ddf = ddf)
      
      
      
      result$ttest <- list(df=tsummary[, "df"], tvalue=tsummary[, "t value"], 
                           tpvalue=tsummary[, "p-value"])
      return(result)
    }
    
    
    #calculate lsmeans of differences of LSMEANS of the final model
    if(isLSMEANS || isDiffLSMEANS)
    {
      if(isLSMEANS)
      {
        lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, 
                                   test.effs=test.effs,
                                   lsmeansORdiff=TRUE, 
                                   l.lmerTest.private.contrast)
        result$lsmeans.table <- lsmeans.tab$summ.data
        result$diffs.lsmeans.table <- NULL
      }
      if(isDiffLSMEANS)
      {
        lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, 
                                   test.effs=test.effs,
                                   lsmeansORdiff=FALSE, 
                                   l.lmerTest.private.contrast)
        result$diffs.lsmeans.table <- lsmeans.tab$summ.data
        result$lsmeans.table <- NULL
      }
      return(result)
    }
    
    
    # Calculate  F-test with Satterthwaite's approximation
    # create X design matrix for fixed effects
    X.design.list <- createDesignMat(model,data)
    X.design <- X.design.list$X.design
    names.design.withLevels <- X.design.list$names.design.withLevels
    
    
   
    ###new code with X.design matrix
    fullCoefs <- rep(0, ncol(X.design))
    fullCoefs <- setNames(fullCoefs, names.design.withLevels) 
    if("(Intercept)" %in% names.design.withLevels)
      names(fullCoefs)[1] <- "(Intercept)"
    fullCoefs[names(rho$fixEffs)] <- rho$fixEffs
    rho$nums.Coefs <- which(names(fullCoefs) %in% names(rho$fixEffs))
   
    rho$nums.Coefs <- setNames(rho$nums.Coefs, names(fullCoefs[rho$nums.Coefs]))
    #define the terms that are to be tested
    test.terms <- attr(terms(model),"term.labels")
    
    #initialize anova table
    if(is.first.anova)
    {
      anova.table <- initAnovaTable(model, test.terms, reduce.fixed)
      is.first.anova <- FALSE
      elim.num <- 1
    }
    
    # calculate general set matrix for type 3 hypothesis
    if(type == 3)
      L <- calcGeneralSetForHypothesis(X.design, rho)  
    
    
    # calculate type 1 hypothesis matrices for each term   
    if(type == 1 || type == 2)
    {
      X <- X.design
      p <- ncol(X)
      XtX <- crossprod(X)
      U <- doolittle(XtX)$U
      d <- diag(U)
      for(i in 1:nrow(U))
        if(d[i] > 0) U[i, ] <- U[i, ] / d[i]
      L <- U
    }
  
      resultFpvalueSS <- llply(test.terms, calcFpvalueMAIN, L = L, 
                               X.design = X.design,
                               fullCoefs = fullCoefs, model = model, rho = rho, 
                               ddf = ddf, type = type)

       
    #fill anova table
    anova.table <- fillAnovaTable(resultFpvalueSS,  anova.table)
    
   
    if(!reduce.fixed)             
      break     
    else
    {
      resNSelim <- elimNSFixedTerm(model, anova.table, data, alpha.fixed, elim.num,
                                   l.lmerTest.private.contrast, keep.effs$fixedeffs)
      if(is.null(resNSelim))
        break
      else
      {
        model <- resNSelim$model
        mf.final <- update.formula(formula(model),formula(model))
        model <- updateModel(model, mf.final, getREML(model), 
                             l.lmerTest.private.contrast)        
        anova.table <- updateAnovaTable(resNSelim)
        elim.num <- elim.num+1
       
      }        
    }      
    
  }
  
  #convert anova table to data frame
  anova.table <- as.data.frame(anova.table)
  anova.table$NumDF <- as.integer(anova.table$NumDF)
  if(!is.null(anova.table$elim.num))
    anova.table$elim.num <- as.integer(anova.table$elim.num)
  
  if(isTotal || isAnova)
  {
    result$anova.table <- anova.table
    if(isAnova)
      return(result)
  }
  
  #if in step function least squares means of diffs of LSMEANS are required
  if(lsmeans.calc)
  {
    lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, 
                               test.effs = test.effs,
                               lsmeansORdiff = TRUE, 
                               l.lmerTest.private.contrast)
    result$lsmeans.table <- lsmeans.tab$summ.data
  }
  else
  {
    result$lsmeans.table <- NULL
  }
  if(difflsmeans.calc)
  {
    lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, 
                               test.effs = test.effs, 
                               lsmeansORdiff=FALSE, 
                               l.lmerTest.private.contrast)
    result$diffs.lsmeans.table <- lsmeans.tab$summ.data
  }
  else
  {
    result$diffs.lsmeans.table <- NULL
  }
  
    
  #format anova.table and random.table according to elim.num column
  result$anova.table <- formatElimNumTable(result$anova.table) 
  result$rand.table <- formatElimNumTable(result$rand.table) 
  
  #update final model
  mf.final <- update.formula(formula(model),formula(model))
  model <- updateModel(model, mf.final, getREML(model), contr)
  
  #save model
  if(inherits(model, "merMod"))
    model <- as(model,"merModLmerTest")
  
  result$model <- model
  return(result)
}



