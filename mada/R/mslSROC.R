## MosesShapiroLittenberg
mslSROC <- function(data = NULL, subset=NULL,
  TP="TP", FN="FN", FP="FP", TN="TN", 
  fpr = NULL, extrapolate = FALSE, 
  correction = 0.5, correction.control = "all",
  add = FALSE, lty = 1, lwd = 1, col = 1, ...){
  
  alphasens <- 1
  alphafpr <- 1
  
  stopifnot(is.numeric(correction), 0 <= correction,  
            correction.control %in% c("all", "single", "none"),
            0 <= alphasens, alphasens <= 2, 0 <= alphafpr, alphafpr <= 2,
            is.numeric(TP) | (is.character(TP) & length(TP) == 1),
            is.numeric(FP) | (is.character(FP) & length(FP) == 1),
            is.numeric(TN) | (is.character(TN) & length(TN) == 1),
            is.numeric(FN) | (is.character(FN) & length(FN) == 1),
            is.logical(extrapolate))
  
  
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
  SENS<-TP/number.of.pos
  FPR <- FP/number.of.neg
  
  if(is.null(fpr)){
    fpr <- 1:99/100
    if(extrapolate){bound = c(0,1)}
    if(!extrapolate){bound = c(min(FPR), max(FPR))}
    fpr <- fpr[cut(fpr,bound, "withinbound") == "withinbound"]  
  }  
    
  senstrafo <- function(x){return(talpha(alphasens)$linkfun(x))}
  fprtrafo <- function(x){return(talpha(alphafpr)$linkfun(x))}
  sensinvtrafo <- function(x){return(talpha(alphasens)$linkinv(x))}
  
  z <- fprtrafo(FPR)
  y <- senstrafo(SENS)  
  D <- y - z
  S <- y + z
  fit <- lm(D~S) 
  A1 <- fit$coefficients[1]
  B1 <- fit$coefficients[2]
  A2 <-  A1/(1-B1)
  B2  <- (1+B1)/(1-B1)

  sens <- sensinvtrafo(A2+B2*fprtrafo(fpr))
  
  if(!add){plot(fpr, sens, type = "l", 
                lty = lty, lwd = lwd, col = col, ...)
           }else{
             lines(fpr, sens, 
                   lty = lty, lwd = lwd, col = col, ...)
           }
  
  return(invisible(list(fpr = fpr, sens = sens, A1 = A1, B1 = B1,
                        A2 = A2, B2 = B2)))
}


