# auxiliar function: estimate + CI for PPV and NPV when casecontrol=TRUE
mercaldo <- function(tab, p, which, conf.level=0.95){
  TP <- tab[1,1]
  FP <- tab[1,2]
  FN <- tab[2,1]
  TN <- tab[2,2]
  
  # Sensitivity
  Se <- TP/(TP + FN)
  # Specificity
  Sp <- TN/(FP + TN)
  
  z <- qnorm(1-(1-conf.level)/2)
  if (which=="PPV"){
    # Positive predictive value
    PPV <- Se*p/(Se*p+(1-Sp)*(1-p))
    # CI for PPV
    logit.PPV <- log(p*Se/((1-Sp)*(1-p)))
    var.logit.PPV <- ((1-Se)/Se)*(1/(TP+FN))+(Sp/(1-Sp))*(1/(FP+TN))
    PPV.ll <- (exp(logit.PPV - z*sqrt(var.logit.PPV)))/(1 + exp(logit.PPV - z*sqrt(var.logit.PPV)))
    PPV.ul <- (exp(logit.PPV + z*sqrt(var.logit.PPV)))/(1 + exp(logit.PPV + z*sqrt(var.logit.PPV)))
    PPV.ic <- c(PPV, PPV.ll, PPV.ul)
    return(PPV.ic)
  }else if (which=="NPV"){
    # Negative predictive value
    NPV <- Sp*(1-p)/(Sp*(1-p)+(1-Se)*p)
    # CI for NPV
    logit.NPV <- log((1-p)*Sp/((1-Se)*p))
    var.logit.NPV <- (Se/(1-Se))*(1/(TP+FN))+((1-Sp)/Sp)*(1/(FP+TN))
    NPV.ll <- (exp(logit.NPV - z*sqrt(var.logit.NPV)))/(1 + exp(logit.NPV - z*sqrt(var.logit.NPV)))
    NPV.ul <- (exp(logit.NPV + z*sqrt(var.logit.NPV)))/(1 + exp(logit.NPV + z*sqrt(var.logit.NPV)))
    NPV.ic <- c(NPV, NPV.ll, NPV.ul)
    return(NPV.ic)
  }else{
    stop("'which' must be 'PPV' or 'NPV'")
  }  
}

# auxiliar function: estimate + CI for OR
OR <- function(tab, conf.level=0.95){
  TP <- tab[1,1]
  FP <- tab[1,2]
  FN <- tab[2,1]
  TN <- tab[2,2]
  
  # Sensitivity
  Se <- TP/(TP + FN)
  # Specificity
  Sp <- TN/(FP + TN)
  
  z <- qnorm(1-(1-conf.level)/2)
  OR <- (Se/(1-Se))/((1-Sp)/Sp)
  var.ln.OR <- 1/TP+1/FP+1/FN+1/TN
  OR.ic <- c(OR, OR*exp(-z*sqrt(var.ln.OR)), OR*exp(z*sqrt(var.ln.OR)))
  return(OR.ic)
}

# auxiliar function: estimate + CI for LR+ and LR-
LR <- function(tab, which, conf.level=0.95){
  TP <- tab[1,1]
  FP <- tab[1,2]
  FN <- tab[2,1]
  TN <- tab[2,2]
  
  # Sensitivity
  Se <- TP/(TP + FN)
  # Specificity
  Sp <- TN/(FP + TN)
  
  z <- qnorm(1-(1-conf.level)/2)
  if (which=="+"){
    LRP <- Se/(1-Sp)
    var.ln.LRP <- (1-Se)/TP+Sp/FP
    LRP.ic <- c(LRP, LRP*exp(-z*sqrt(var.ln.LRP)), LRP*exp(z*sqrt(var.ln.LRP)))
    return(LRP.ic)
  }else if (which=="-"){
    LRN <- (1-Se)/Sp
    var.ln.LRN <- Se/FN+(1-Sp)/TN
    LRN.ic <- c(LRN, LRN*exp(-z*sqrt(var.ln.LRN)), LRN*exp(z*sqrt(var.ln.LRN)))
    return(LRN.ic)
  }else{
    stop("'which' must be '+' or '-'")
  }
}



# diag: diagnostic test
diagnostic <- function(tab, method=c("par", "exact"), casecontrol=FALSE, p=NULL, conf.level=0.95){
  
  # arguments:
  #  tab: 2x2-dimensional table in the following form:
  #         TP FP
  #         FN TN
  #  method: method for computing the CIs for Se, Sp, PPV, NPV, accuracy and error rate.
  #          The user can choose between "par" (parametric) and "exact" (exact).
  #          Default, "par".
  #  casecontrol: Were data collected in a case-control study?
  #               Default, FALSE.
  #  p: disease prevalence (only when casecontrol=T).
  #  conf.level: level confidence for the CIs. Default, 95%.
  
  # entry control
  method <- match.arg(method)
  if (!is.table(tab) & !is.matrix(tab)){
    stop("'tab' must be a table or a matrix")
  }
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                                 conf.level < 0 || conf.level > 1)){
    stop("'conf.level' must be a single number between 0 and 1")
  }
  if (!all(dim(tab)==2)){
    stop("'tab' must be a 2x2 table or a matrix")
  }
  
  # values in tab
  TP <- tab[1,1]
  FP <- tab[1,2]
  FN <- tab[2,1]
  TN <- tab[2,2]
  
  # Sensitivity (Se)
  Se <- TP/(TP + FN)
  # Specificity (Sp)
  Sp <- TN/(FP + TN)
  
  # CI for Se and Sp
  if (method=="par"){
    Se.ic <- c(Se, prop.test(TP, TP + FN, conf.level=conf.level)$conf.int)
    Sp.ic <- c(Sp, prop.test(TN, FP + TN, conf.level=conf.level)$conf.int)
  }else if (method=="exact"){
    Se.ic <- c(Se, binom.test(TP, TP + FN, conf.level=conf.level)$conf.int)
    Sp.ic <- c(Sp, binom.test(TN, FP + TN, conf.level=conf.level)$conf.int)
  }
  
  # PPV and NPV
  if (!casecontrol){
    # Positive predictive value
    PPV <- TP/(TP + FP)
    # Negative predictive value
    NPV <- TN/(TN + FN)
    # CI for PPV and NPV
    if (method=="par"){
      PPV.ic <- c(PPV, prop.test(TP, TP + FP, conf.level=conf.level)$conf.int)
      NPV.ic <- c(NPV, prop.test(TN, TN + FN, conf.level=conf.level)$conf.int)
    }else if (method=="exact"){
      PPV.ic <- c(PPV, binom.test(TP, TP + FP, conf.level=conf.level)$conf.int)
      NPV.ic <- c(NPV, binom.test(TN, TN + FN, conf.level=conf.level)$conf.int)
    }
  }else{
    if (is.null(p)){
      stop("prevalence 'p' must be provided when casecontrol=TRUE")
    }else if (p<=0 | p>=1){
      stop("'p' must be in (0,1)")
    }else{
      PPV.ic <- mercaldo(tab, p, "PPV", conf.level)
      NPV.ic <- mercaldo(tab, p, "NPV", conf.level)
      # applying continuity correction when PPV or NPV are estimated as 1 
      if (PPV.ic[1]==1){
        PPV.ic <- mercaldo(tab + 0.5, p, "PPV", conf.level)
      }
      if (NPV.ic[1]==1){
        NPV.ic <- mercaldo(tab + 0.5, p, "NPV", conf.level)
      }
    }
  }
  
  # LR+
  if (Sp==1 | Se==0){
    # applying continuity correction
    LRP.ic <- LR(tab+0.5, "+", conf.level)
  }else{
    LRP.ic <- LR(tab, "+", conf.level)
  }  
  
  # LR-
  if (Sp==0 | Se==1){
    # applying continuity correction
    LRN.ic <- LR(tab+0.5, "-", conf.level)
  }else{
    LRN.ic <- LR(tab, "-", conf.level)
  }
  
  # Odds ratio
  if (any(tab==0)){
    # applying continuity correction
    OR.ic <- OR(tab+0.5, conf.level)
  }else{
    OR.ic <- OR(tab, conf.level)
  }
  
  # Youden index
  YI <- Se + Sp - 1 
  var.YI <- Se*(1-Se)/(TP+FN) + Sp*(1-Sp)/(FP+TN)
  z <- qnorm(1-(1-conf.level)/2)
  YI.ic <- c(YI, YI - z*sqrt(var.YI), YI + z*sqrt(var.YI))
  
  # Accuracy (probability of a correct test result)
  Acc <- (TP+TN)/(TP+FP+FN+TN)
  # CI for accuracy
  if (method=="par"){
    Acc.ic <- c(Acc, prop.test(TP+TN, TP+FP+FN+TN, conf.level=conf.level)$conf.int)
  }else if (method=="exact"){
    Acc.ic <- c(Acc, binom.test(TP+TN, TP+FP+FN+TN, conf.level=conf.level)$conf.int)
  }
  
  # Error rate
  ER <- (FP + FN)/(TP+TN+FP+FN)
  # CI for error rate
  if (method=="par"){
    ER.ic <- c(ER, prop.test(FP+FN, TP+FP+FN+TN, conf.level=conf.level)$conf.int)
  }else if (method=="exact"){
    ER.ic <- c(ER, binom.test(FP+FN, TP+FP+FN+TN, conf.level=conf.level)$conf.int)
  }
  
  # table of results
  res <- rbind(Se.ic, Sp.ic, PPV.ic, NPV.ic, LRP.ic, LRN.ic, OR.ic, YI.ic, Acc.ic,
               ER.ic)
  
  colnames(res) <- c("Estim.", paste("Low.lim(", conf.level*100, "%)", sep=""),
                     paste("Up.lim(", conf.level*100, "%)", sep=""))
  rownames(res) <- c("Sensitivity", "Specificity", "Pos.Pred.Val.", "Neg.Pred.Val.",
                     "LR+", "LR-", "Odds ratio", "Youden index", "Accuracy",
                     "Error rate")
  return(res)
}