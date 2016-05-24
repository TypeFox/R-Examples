getSummary.mclogit <- function(obj,
                               alpha=.05,
                               rearrange=NULL,
                               #as.columns=NULL,
                               ...){
  
  smry <- summary(obj)
  N <- obj$N
  coef <- smry$coefficients
  varPar <- smry$varPar
  
  lower.cf <- qnorm(p=alpha/2,mean=coef[,1],sd=coef[,2])
  upper.cf <- qnorm(p=1-alpha/2,mean=coef[,1],sd=coef[,2])
  coef <- cbind(coef,lower.cf,upper.cf)
  colnames(coef) <- c("est","se","stat","p","lwr","upr")
  if(length(varPar)){
    se.log.varPar <- varPar[,1]*varPar[,2]
    lower.log.varPar <- qnorm(p=alpha/2,mean=log(varPar[,1]),sd=se.log.varPar[2])
    upper.log.varPar <- qnorm(p=1-alpha/2,mean=log(varPar[,1]),sd=se.log.varPar[2])
    varPar <- cbind(varPar,exp(lower.log.varPar),exp(upper.log.varPar))
    colnames(varPar) <- c("est","se","stat","p","lwr","upr")
    rownames(varPar) <- paste("Var(",rownames(varPar),")",sep="")
  }
  if(length(rearrange)){
    coef.grps <- lapply(rearrange,function(ii){
      if(is.character(ii) && !all(ii %in% rownames(coef)))
        stop("coefficient(s) ",dQuote(unname(ii[!(ii %in% rownames(coef))]))," do not exist")
      structure(coef[ii,],
                dimnames=list(names(ii),dimnames(coef)[[2]])
      )
    })
    grp.titles <- names(rearrange)
    coef.grps <- do.call(memisc::collect,coef.grps)
    coef <- array(NA,dim=c(
      dim(coef.grps)[1] + NROW(varPar),
      dim(coef.grps)[2],
      dim(coef.grps)[3]
    ))
    coef[seq(dim(coef.grps)[1]),,] <- coef.grps
    if(length(varPar))
      coef[dim(coef.grps)[1]+seq(nrow(varPar)),,1] <- varPar
    dimnames(coef) <- list(
      c(dimnames(coef.grps)[[1]],rownames(varPar)),
      dimnames(coef.grps)[[2]],
      grp.titles
    )
  }
  else {
    .coef <- coef
    coef <- matrix(NA,nrow=nrow(.coef)+NROW(varPar),ncol=ncol(.coef))
    coef[seq(nrow(.coef)),] <- .coef
    if(length(varPar))
      coef[nrow(.coef)+seq(nrow(varPar)),] <- varPar
    rownames(coef) <- c(rownames(.coef),rownames(varPar))
    colnames(coef) <- colnames(.coef)
  }
  
  
  phi <- smry$phi
  LR <- smry$null.deviance - smry$deviance
  df <- obj$model.df
  deviance <- deviance(obj)
  
  
  if(df > 0){
    p <- pchisq(LR,df,lower.tail=FALSE)
    L0.pwr <- exp(-smry$null.deviance/N)
    LM.pwr <- exp(-smry$deviance/N)
    
    McFadden <- 1- smry$deviance/smry$null.deviance
    Cox.Snell <- 1 - exp(-LR/N)
    Nagelkerke <- Cox.Snell/(1-L0.pwr)
  }
  else {
    LR <- NA
    df <- NA
    p <- NA
    McFadden <- NA
    Cox.Snell <- NA
    Nagelkerke <- NA
  }
  
  ll <- obj$ll
  AIC <- AIC(obj)
  BIC <- AIC(obj,k=log(N))
  sumstat <- c(
    phi         = phi,
    LR             = LR,
    df         = df,
    #p             = p,
    logLik        = ll,
    deviance      = deviance,
    McFadden      = McFadden,
    Cox.Snell       = Cox.Snell,
    Nagelkerke    = Nagelkerke,
    AIC           = AIC,
    BIC           = BIC,
    N             = N
  )
  
  #coef <- apply(coef,1,applyTemplate,template=coef.template)
  
  #sumstat <- drop(applyTemplate(sumstat,template=sumstat.template))
  list(coef=coef,sumstat=sumstat)
}
