



estimate.ATE <- function(pscore.formula, pscore.family,
                         outcome.formula.t, outcome.formula.c,
                         outcome.family, treatment.var,
                         data=NULL,
                         divby0.action=c("fail", "truncate", "discard"),
                         divby0.tol=1e-8,
                         nboot=501, variance.smooth.deg=1,
                         variance.smooth.span=0.75,
                         var.gam.plot=TRUE, suppress.warnings=TRUE, ...){

  
  if (is.null(data)){
    stop("'data' must be specified in call to 'estimate.ATE'\n")
  }

  if (nboot < 0){
    nboot <- 0
  }
  
  call <- match.call()
  arg.list <- as.list(call)
  
  if (("zzzzzbsfitzzzzz" %in% names(arg.list))){
    zzzzzbsfitzzzzz <- TRUE ## don't calculate asymp var
  }
  else{
    zzzzzbsfitzzzzz <- FALSE ## calculate asymp var
  }
  if (variance.smooth.span <=0 | variance.smooth.deg < 0){
    zzzzzbsfitzzzzz <- TRUE  ## don't calculate asymp var
  }

  
  divby0.action <- match.arg(divby0.action)

  if (max(is.na(data)) > 0){
    stop("'data' contains NAs. Remove NAs and call estimate.ATE again.\n")
  }



  
  #pscore.model.frame <- model.frame(pscore.formula, data=data, ...)
  #treatment.vec <- model.response(pscore.model.frame)
  #treatment.var <- colnames(pscore.model.frame)[1]
 

  treatment.vec <- data[, treatment.var]
  
  if (is.factor(treatment.vec)){
    treatment.values <- levels(treatment.vec)
  }
  else{
    treatment.values <- sort(unique(treatment.vec))
  }



  
  treated.data <- data[treatment.vec==treatment.values[2],]
  control.data <- data[treatment.vec==treatment.values[1],]

  
  
  n.treated.pre <- nrow(treated.data)
  n.control.pre <- nrow(control.data)

  
  if (suppress.warnings){
    gam.ps <- suppressWarnings(gam(pscore.formula, family=pscore.family,
                  data=data, na.action="na.fail", ...))
  }
  else{
    gam.ps <- gam(pscore.formula, family=pscore.family,
                  data=data, na.action="na.fail", ...)
  }


  
  pscore.probs <- predict(gam.ps, newdata=data, type="response", ...)


  

  
  pscores.pre <- pscore.probs
    
  truncated.indic <- rep(FALSE, nrow(data))
  discarded.indic <- rep(FALSE, nrow(data))
  n.treated.post <- n.treated.pre
  n.control.post <- n.control.pre
  if(min(pscore.probs) <= divby0.tol || max(pscore.probs) >= (1-divby0.tol)){
    if (divby0.action == "fail"){
      stop("\nCannot compute AIPW estimate because some\nprobabilities of treatment are numerically 0 and/or 1\n\n")
    }
    if (divby0.action == "truncate"){
      truncated.indic[pscore.probs <= divby0.tol] <- TRUE
      pscore.probs[pscore.probs <= divby0.tol] <- divby0.tol
      truncated.indic[pscore.probs >= (1-divby0.tol)] <- TRUE
      pscore.probs[pscore.probs >= (1-divby0.tol)] <- (1-divby0.tol)
      n.treated.post <- sum((treatment.vec==treatment.values[2])[!truncated.indic])
      n.control.post <- sum((treatment.vec==treatment.values[1])[!truncated.indic])
    }
    if (divby0.action == "discard"){
      discarded.indic[pscore.probs <= divby0.tol] <- TRUE
      discarded.indic[pscore.probs >= (1-divby0.tol)] <- TRUE
      pscore.probs <- pscore.probs[!discarded.indic]
      treatment.vec <- treatment.vec[!discarded.indic]
      data <- data[!discarded.indic,]
      treated.data <- data[treatment.vec==treatment.values[2],]
      control.data <- data[treatment.vec==treatment.values[1],]
      n.treated.post <- nrow(treated.data)
      n.control.post <- nrow(control.data)
    }    
  }




  outcome.vec <- model.response(model.frame(outcome.formula.t, data=data))


  
  if (suppress.warnings){
    gam.t <- suppressWarnings(gam(outcome.formula.t, family=outcome.family,
                 data=treated.data, na.action="na.fail", ...))
    gam.c <- suppressWarnings(gam(outcome.formula.c, family=outcome.family,
                 data=control.data, na.action="na.fail", ...))
    outcome.treated.expectation <- suppressWarnings(predict(gam.t, newdata=data,
                                           type="response", ...))  
    outcome.control.expectation <- suppressWarnings(predict(gam.c, newdata=data,
                                           type="response", ...))    
    res.c <- suppressWarnings(residuals(gam.c, type="response"))
    res.t <- suppressWarnings(residuals(gam.t, type="response"))
  }
  else{
    gam.t <- gam(outcome.formula.t, family=outcome.family,
                 data=treated.data, na.action="na.fail", ...)
    gam.c <- gam(outcome.formula.c, family=outcome.family,
                 data=control.data, na.action="na.fail", ...)
    outcome.treated.expectation <- predict(gam.t, newdata=data,
                                           type="response", ...)  
    outcome.control.expectation <- predict(gam.c, newdata=data,
                                           type="response", ...)
    res.c <- residuals(gam.c, type="response")
    res.t <- residuals(gam.t, type="response")
  }
  
  

  ## number of observations
  n <- length(treatment.vec)


  term1 <- rep(0, n)
  term2 <- rep(0, n)
  term3 <- rep(0, n)
  term4 <- rep(0, n)
  control.indic <- treatment.vec == treatment.values[1]
  treated.indic <- treatment.vec == treatment.values[2]


  term1 <- outcome.vec/pscore.probs
  term1[control.indic] <- 0

  term2 <- ((1 - pscore.probs) * outcome.treated.expectation) /
    pscore.probs
  term2[control.indic] <- ((0 - pscore.probs[control.indic]) * outcome.treated.expectation[control.indic]) /pscore.probs[control.indic]

  term3 <-  outcome.vec / (1 - pscore.probs)
  term3[treated.indic] <- 0

  term4 <- ((1 - pscore.probs) * outcome.control.expectation) / (1 - pscore.probs)
  term4[control.indic] <- ((0 - pscore.probs[control.indic]) * outcome.control.expectation[control.indic]) / (1 - pscore.probs[control.indic])

  
##  for (i in 1:n){
##    if (treatment.vec[i] != treatment.values[1]){ ## treated
##      term1[i] <-  outcome.vec[i]/pscore.probs[i]
##      term2[i] <- ((1 - pscore.probs[i]) * outcome.treated.expectation[i]) /
##        pscore.probs[i] 
##      ## term3[i] <-  0
##      term4[i] <- ((1 - pscore.probs[i]) * outcome.control.expectation[i]) /
##        (1 - pscore.probs[i])
##    }
##    else{ ## control
##      ##term1[i] <- 0
##      term2[i] <- ((0 - pscore.probs[i]) * outcome.treated.expectation[i]) /
##        pscore.probs[i] 
##      term3[i] <- outcome.vec[i] / (1 - pscore.probs[i])
##      term4[i] <- ((0 - pscore.probs[i]) * outcome.control.expectation[i]) /
##        (1 - pscore.probs[i])       
##    }    
##  } ## end i in 1:n loop

  ATE.AIPW.hat <- (1/n) * (sum(term1) - sum(term2) - sum(term3) - sum(term4))

    
  I.hat <- term1 - term2 - term3 - term4 - ATE.AIPW.hat
  ATE.AIPW.sand.var <- (1 / n^2) * sum(I.hat^2)
  ATE.AIPW.sand.se <- sqrt(ATE.AIPW.sand.var)





  ## regression
  ATE.reg.hat <- mean(outcome.treated.expectation) -
    mean(outcome.control.expectation) 
  


  ## IPW
  weight.norm.1 <- 1/sum(treated.indic/pscore.probs)
  weight.norm.3 <- 1/sum((1-treated.indic)/(1-pscore.probs))
  ATE.IPW.hat <- weight.norm.1 * sum(term1) - weight.norm.3 * sum(term3)



  ## Calculate estimated large sample SEs
  res.data.c <- data.frame(res=res.c,
                           pscore.probs=pscore.probs[control.indic])
  res.data.t <- data.frame(res=res.t,
                           pscore.probs=pscore.probs[treated.indic])



  if (length(unique(pscore.probs)) > 3){
    var.formula <- as.formula(substitute((res^2) ~ lo(pscore.probs,
                                                      degree=variance.smooth.deg,
                                                      span=variance.smooth.span),
                                         list(variance.smooth.deg=variance.smooth.deg, variance.smooth.span=variance.smooth.span)))
  }
  else{
    var.formula <- as.formula(substitute((res^2) ~ pscore.probs,
                                         list(variance.smooth.deg=variance.smooth.deg, variance.smooth.span=variance.smooth.span)))
  }
  
#  var.formula <- as.formula((res^2) ~ lo(pscore.probs,
#                                         degree=2, span=.75))
  
  if(!zzzzzbsfitzzzzz){
    if (suppress.warnings){
      var.gam.c <- suppressWarnings(gam(var.formula, data=res.data.c,
                                       family=gaussian(log)))
      var.gam.t <- suppressWarnings(gam(var.formula, data=res.data.t,
                                        family=gaussian(log)))
    }
    else{
      var.gam.c <- gam(var.formula, data=res.data.c, family=gaussian(log))
      var.gam.t <- gam(var.formula, data=res.data.t, family=gaussian(log))
    }
      
    if(var.gam.plot){
      par(mfrow=c(2,1))
      plot(res.data.c$pscore.probs, res.data.c$res^2,
           xlab="Propensity Score", ylab="Squared Residuals",
           main="Conditional Variance under Control")
      ord.c <- order(res.data.c$pscore.probs)
      lines(res.data.c$pscore.probs[ord.c],
            predict(var.gam.c, type="response")[ord.c], col="red", lwd=2)
      
      
      plot(res.data.t$pscore.probs, res.data.t$res^2,
           xlab="Propensity Score", ylab="Squared Residuals",
           main="Conditional Variance under Treatment")
      ord.t <- order(res.data.t$pscore.probs)
      lines(res.data.t$pscore.probs[ord.t],
            predict(var.gam.t, type="response")[ord.t], col="red", lwd=2)
      
    }


    cond.var.c <- suppressWarnings(predict(var.gam.c,
                                           newdata=data, type="response"))
    cond.var.t <- suppressWarnings(predict(var.gam.t,
                                           newdata=data, type="response"))
  }
  else{
    cond.var.c <- rep(NA, n)
    cond.var.t <- rep(NA, n)
  }
  
  ATE.AIPW.asymp.var <- (1/(n^2))*sum(
                                      cond.var.t/pscore.probs +
                                      cond.var.c/(1 - pscore.probs) +
                                      (outcome.treated.expectation -
                                       outcome.control.expectation -
                                       ATE.AIPW.hat)^2
                                      )
  
  ATE.reg.asymp.var <- (1/(n^2))*sum(
                                     cond.var.t/pscore.probs +
                                     cond.var.c/(1 - pscore.probs) +
                                     (outcome.treated.expectation -
                                      outcome.control.expectation -
                                      ATE.reg.hat)^2
                                     )
  
  ATE.IPW.asymp.var <- (1/(n^2))*sum(
                                     cond.var.t/pscore.probs +
                                     cond.var.c/(1 - pscore.probs) +
                                     (outcome.treated.expectation -
                                      outcome.control.expectation -
                                      ATE.IPW.hat)^2
                                     )
  
  
  ATE.AIPW.asymp.se <- sqrt(ATE.AIPW.asymp.var)
  ATE.reg.asymp.se <- sqrt(ATE.reg.asymp.var)
  ATE.IPW.asymp.se <- sqrt(ATE.IPW.asymp.var)
  



  ## calculate bootstrap SEs
  ATE.AIPW.bs <- rep(NA, nboot)
  ATE.reg.bs <- rep(NA, nboot)
  ATE.IPW.bs <- rep(NA, nboot)
  ATE.AIPW.bs.se <- NA
  ATE.reg.bs.se <- NA
  ATE.IPW.bs.se <- NA
  if (nboot > 0){
    
    for (biter in 1:nboot){
      boot.inds <- sample(1:n, n, replace=TRUE)
      data.bs <- data[boot.inds,]
      
      out <- estimate.ATE(pscore.formula=pscore.formula,
                          pscore.family=pscore.family,
                          outcome.formula.t=outcome.formula.t,
                          outcome.formula.c=outcome.formula.c,
                          outcome.family=outcome.family,
                          treatment.var=treatment.var,
                          variance.smooth.deg=variance.smooth.deg,
                          variance.smooth.span=variance.smooth.span,
                          data=data.bs,
                          divby0.action=divby0.action,
                          divby0.tol=divby0.tol, nboot=0, var.gam.plot=FALSE,
                          suppress.warnings=suppress.warnings,
                          zzzzzbsfitzzzzz=TRUE)
      
      ATE.AIPW.bs[biter] <- out$ATE.AIPW.hat
      ATE.reg.bs[biter] <- out$ATE.reg.hat
      ATE.IPW.bs[biter] <- out$ATE.IPW.hat
    }

    
    ATE.AIPW.bs.se <- sd(ATE.AIPW.bs)
    ATE.reg.bs.se <- sd(ATE.reg.bs)
    ATE.IPW.bs.se <- sd(ATE.IPW.bs)
  }
  

  
  
  return(structure(list(ATE.AIPW.hat=ATE.AIPW.hat,
                        ATE.reg.hat=ATE.reg.hat,
                        ATE.IPW.hat=ATE.IPW.hat,
                        ATE.AIPW.sand.SE=ATE.AIPW.sand.se,
                        ATE.AIPW.asymp.SE=ATE.AIPW.asymp.se,
                        ATE.reg.asymp.SE=ATE.reg.asymp.se,
                        ATE.IPW.asymp.SE=ATE.IPW.asymp.se,
                        ATE.AIPW.bs.SE=ATE.AIPW.bs.se,
                        ATE.reg.bs.SE=ATE.reg.bs.se,
                        ATE.IPW.bs.SE=ATE.IPW.bs.se,
                        ATE.AIPW.bs=ATE.AIPW.bs,
                        ATE.reg.bs=ATE.reg.bs,
                        ATE.IPW.bs=ATE.IPW.bs,
                        gam.t=gam.t,
                        gam.c=gam.c,
                        gam.ps=gam.ps,
                        truncated.indic=truncated.indic,
                        discarded.indic=discarded.indic,
                        treated.value=treatment.values[2],
                        control.value=treatment.values[1],
                        treatment.var=treatment.var,
                        n.treated.prediscard=n.treated.pre,
                        n.control.prediscard=n.control.pre,
                        n.treated.postdiscard=n.treated.post,
                        n.control.postdiscard=n.control.post,
                        pscores.prediscard=pscores.pre,
                        pscores.postdiscard=pscore.probs,
                        cond.var.t=cond.var.t,
                        cond.var.c=cond.var.c,
                        call=call,
                        data=data),
                   class="CausalGAM"
                   ))
  
}



"print.CausalGAM" <- function(x, ...){
  cat("###################################################################\n")
  cat("AIPW Estimator:\n")
  cat("-------------------------------------------------------------------\n")
  cat("Estimated ATE: ", round(x$ATE.AIPW.hat, 4), "\n")
  cat("                     SE            z-statistic            Pr(>|z|)\n")
  pval <- 2*(1-pnorm(abs(x$ATE.AIPW.hat/x$ATE.AIPW.sand.SE)))
  if (pval >= 0.0001){
    cat("Emp. Sandwich     ", round(x$ATE.AIPW.sand.SE, 4), "          ",
        round(x$ATE.AIPW.hat/x$ATE.AIPW.sand.SE, 4), "             ",
        round(pval, 4), "\n")
  }
  else{
    cat("Emp. Sandwich     ", round(x$ATE.AIPW.sand.SE, 4), "          ",
        round(x$ATE.AIPW.hat/x$ATE.AIPW.sand.SE, 4), "             ",
        "< 0.0001", "\n")    
  }
  if (!is.na(x$ATE.AIPW.asymp.SE)){
    pval <- 2*(1-pnorm(abs(x$ATE.AIPW.hat/x$ATE.AIPW.asymp.SE)))
    if (pval >= 0.0001){
      cat("Estim. Asymptotic ", round(x$ATE.AIPW.asymp.SE, 4), "          ",
          round(x$ATE.AIPW.hat/x$ATE.AIPW.asymp.SE, 4), "             ",
          round(pval, 4), "\n")
    }
    else{
      cat("Estim. Asymptotic ", round(x$ATE.AIPW.asymp.SE, 4), "          ",
          round(x$ATE.AIPW.hat/x$ATE.AIPW.asymp.SE, 4), "             ",
          "< 0.0001", "\n")    
    }
  }
  if (!is.na(x$ATE.AIPW.bs.SE)){
    pval <- 2*(1-pnorm(abs(x$ATE.AIPW.hat/x$ATE.AIPW.bs.SE)))
    if (pval >= 0.0001){
      cat("Bootstrap         ", round(x$ATE.AIPW.bs.SE, 4), "          ",
          round(x$ATE.AIPW.hat/x$ATE.AIPW.bs.SE, 4), "             ",
          round(pval, 4), "\n")
    }
    else{
      cat("Bootstrap         ", round(x$ATE.AIPW.bs.SE, 4), "          ",
          round(x$ATE.AIPW.hat/x$ATE.AIPW.bs.SE, 4), "             ",
          "< 0.0001", "\n")    
    }    
  }
  cat("###################################################################\n\n")

  
  cat("###################################################################\n")
  cat("IPW Estimator:\n")
  cat("-------------------------------------------------------------------\n")
  cat("Estimated ATE: ", round(x$ATE.IPW.hat, 4), "\n")
  cat("                     SE            z-statistic            Pr(>|z|)\n")
  if (!is.na(x$ATE.IPW.asymp.SE)){
    pval <- 2*(1-pnorm(abs(x$ATE.IPW.hat/x$ATE.IPW.asymp.SE)))
    if (pval >= 0.0001){
      cat("Estim. Asymptotic ", round(x$ATE.IPW.asymp.SE, 4), "          ",
          round(x$ATE.IPW.hat/x$ATE.IPW.asymp.SE, 4), "             ",
          round(pval, 4), "\n")
    }
    else{
      cat("Estim. Asymptotic ", round(x$ATE.IPW.asymp.SE, 4), "          ",
          round(x$ATE.IPW.hat/x$ATE.IPW.asymp.SE, 4), "             ",
          "< 0.0001", "\n")    
    }
  }
  if (!is.na(x$ATE.IPW.bs.SE)){
    pval <- 2*(1-pnorm(abs(x$ATE.IPW.hat/x$ATE.IPW.bs.SE)))
    if (pval >= 0.0001){
      cat("Bootstrap         ", round(x$ATE.IPW.bs.SE, 4), "          ",
          round(x$ATE.IPW.hat/x$ATE.IPW.bs.SE, 4), "             ",
          round(pval, 4), "\n")
    }
    else{
      cat("Bootstrap         ", round(x$ATE.IPW.bs.SE, 4), "          ",
          round(x$ATE.IPW.hat/x$ATE.IPW.bs.SE, 4), "             ",
          "< 0.0001", "\n")    
    }    
  }
  cat("###################################################################\n\n")


  cat("###################################################################\n")
  cat("Regression Estimator:\n")
  cat("-------------------------------------------------------------------\n")
  cat("Estimated ATE: ", round(x$ATE.reg.hat, 4), "\n")
  cat("                     SE            z-statistic            Pr(>|z|)\n")
  if (!is.na(x$ATE.reg.asymp.SE)){
    pval <- 2*(1-pnorm(abs(x$ATE.reg.hat/x$ATE.reg.asymp.SE)))
    if (pval >= 0.0001){
      cat("Estim. Asymptotic ", round(x$ATE.reg.asymp.SE, 4), "          ",
          round(x$ATE.reg.hat/x$ATE.reg.asymp.SE, 4), "             ",
          round(pval, 4), "\n")
    }
    else{
      cat("Estim. Asymptotic ", round(x$ATE.reg.asymp.SE, 4), "          ",
          round(x$ATE.reg.hat/x$ATE.reg.asymp.SE, 4), "             ",
          "< 0.0001", "\n")    
    }
  }
  if (!is.na(x$ATE.reg.bs.SE)){
    pval <- 2*(1-pnorm(abs(x$ATE.reg.hat/x$ATE.reg.bs.SE)))
    if (pval >= 0.0001){
      cat("Bootstrap         ", round(x$ATE.reg.bs.SE, 4), "          ",
          round(x$ATE.reg.hat/x$ATE.reg.bs.SE, 4), "             ",
          round(pval, 4), "\n")
    }
    else{
      cat("Bootstrap         ", round(x$ATE.reg.bs.SE, 4), "          ",
          round(x$ATE.reg.hat/x$ATE.reg.bs.SE, 4), "             ",
          "< 0.0001", "\n")    
    }    
  }
  cat("###################################################################\n\n")

  cat("###################################################################\n")
  cat("General Information\n")
  cat("-------------------------------------------------------------------\n")
  cat("Control Value: ", paste(x$treatment.var, "=", x$control.value), "\n")
  cat("Treated Value: ", paste(x$treatment.var, "=", x$treated.value), "\n")
  cat("Number of Discarded Units: ", sum(x$discarded.indic), "\n")
  cat("Number of Truncated Propensity Scores: ", sum(x$truncated.indic), "\n")
  cat("Number of Treated Units Before Discards/Truncations: ", x$n.treated.prediscard, "\n")
  cat("Number of Treated Units not Discarded/Truncated:     ", x$n.treated.postdiscard, "\n")
  cat("Number of Control Units Before Discards/Truncations: ", x$n.control.prediscard, "\n")
  cat("Number of Control Units not Discarded/Truncated:     ", x$n.control.postdiscard, "\n")
  cat("###################################################################\n\n")
  
  
}


