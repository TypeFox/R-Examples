## Ruecker Schumacher
rsSROC <- function(data = NULL, subset=NULL,
  TP="TP", FN="FN", FP="FP", TN="TN", 
  lambda = "from_bivariate",                  
  fpr = NULL, extrapolate = FALSE, plotstudies = FALSE,
  correction = 0.5, correction.control = "all",
  add = FALSE, lty = 1, lwd = 1, col = 1, ...){
  
  
  stopifnot(is.numeric(correction), 0 <= correction,  
            correction.control %in% c("all", "single", "none"),
            is.numeric(TP) | (is.character(TP) & length(TP) == 1),
            is.numeric(FP) | (is.character(FP) & length(FP) == 1),
            is.numeric(TN) | (is.character(TN) & length(TN) == 1),
            is.numeric(FN) | (is.character(FN) & length(FN) == 1),
            is.logical(extrapolate),
            lambda == "from_bivariate" | 
              (is.numeric(lambda) & length(lambda) == 1))
  
  if(!is.null(subset)){
    if(!is.null(data)){data <- data[subset,]}else{
      TP <- TP[subset]
      FP <- FP[subset]
      TN <- TN[subset]
      FN <- FN[subset]
    }
  }
  
  if(!is.null(data) & is.character(c(TP,FP,TN,FN))){
    X <- as.data.frame(data)
    origdata <- data
    TP <- getElement(X,TP)
    FN <- getElement(X,FN)
    FP <- getElement(X,FP)
    TN <- getElement(X,TN)
  }
  
  if(is.null(data) & !is.character(c(TP,FP,TN,FN))){
    origdata <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN)
  }
  
  freqdata <- cbind(TP,FN,FP,TN)
  checkdata(freqdata)
  
  N <- length(TP)  
  
  ## apply continuity correction to _all_ studies if one contains zero
  if(correction.control == "all"){if(any(c(TP,FN,FP,TN) == 0))
  {TP <- TP + correction;
   FN <- FN + correction;
   FP <- FP + correction;
   TN <- TN + correction}}
  if(correction.control == "single"){
    correction = ((((TP == 0)|(FN == 0))|(FP == 0))| (TN == 0)) * 
      correction
    TP <- correction + TP
    FN <- correction + FN
    FP <- correction + FP
    TN <- correction + TN
  }
  
  number.of.pos <- TP + FN
  number.of.neg <- FP + TN
  Se <-TP/number.of.pos
  Sp <- 1 - FP/number.of.neg
  
  if(is.null(fpr)){
    fpr <- 1:99/100
    if(extrapolate){bound = c(0,1)}
    if(!extrapolate){bound = c(min(1-Sp), max(1-Sp))}
    fpr <- fpr[cut(fpr,bound, "withinbound") == "withinbound"]  
  }  
   
  if(lambda == "from_bivariate"){
    fit <- reitsma(data.frame(TP = round(TP), FP = round(FP),
                              FN = round(FN), TN = round(TN)),
                   correction = correction,
                   correction.control = correction.control)  
    Se0 <- expit(fit$coefficients[1,1])  
    Sp0 <- 1-expit(fit$coefficients[1,2]) 
    lambda <- Sp0*(1-Sp0)/(Sp0*(1-Sp0)+Se0*(1-Se0))
  }
  
  y <- logit(Se)
  z <- logit(1-Sp)
  
  bb <- (1-lambda)*Sp*(1-Sp)/Se/(1-Se)/lambda
  aa <- y-bb*z
  
  # The weights are the inverse variances of the log slope estimate
  v.logbb <- 1/(TP) + 1/(FN) + 1/(TN) + 1/(FP)
  w <- 1/v.logbb
  beta2 <- exp(sum(w*log(bb))/sum(w))
  alpha2 <- sum(w*(y-bb*z)/sum(w))
  
  sens <- expit(alpha2 + beta2*log(fpr/(1-fpr)))

  if(!add){plot(fpr, sens, type = "l", 
                lty = lty, lwd = lwd, col = col, ...)
  }else{
    lines(fpr, sens, 
          lty = lty, lwd = lwd, col = col, ...)
  }
  
  if(plotstudies){
    for(i in 1:length(bb)){
      lines(fpr, expit(aa[i] + bb[i]*log(fpr/(1-fpr))))
    }
  }
  
  return(invisible(list(fpr = fpr, sens = sens, lambda = lambda,
                        bb = bb, aa = aa, v.logbb = v.logbb,
                        beta2 = beta2, alpha2 = alpha2)))
}


