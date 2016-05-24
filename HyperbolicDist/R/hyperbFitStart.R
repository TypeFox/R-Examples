hyperbFitStart <- function(x, breaks = NULL,
                           startValues = "BN",
                           ThetaStart = NULL,
                           startMethodSL = "Nelder-Mead",
                           startMethodMoM = "Nelder-Mead", ...)
{
  histData <- hist(x, plot = FALSE, right = FALSE)
  if(is.null(breaks)){
    breaks <- histData$breaks
  }
  midpoints <- histData$mids
  empDens <- ifelse(!is.finite(log(histData$density)), NA, histData$density)
  maxIndex <- order(empDens, na.last = FALSE)[length(empDens)]

  if(length(na.omit(empDens[1:maxIndex]))>1){
    leftAsymptote <- lm(log(empDens)[1:maxIndex]~midpoints[1:maxIndex])$coef
    rightAsymptote <- c(NA,-10*leftAsymptote[2]) # arbitrary large value
  }
  if(length(na.omit(empDens[maxIndex:length(empDens)]))>1){
    rightAsymptote <- lm(log(empDens)[maxIndex:length(empDens)]~
                       midpoints[maxIndex:length(empDens)])$coef
    if(length(na.omit(empDens[1:maxIndex]))<2){
      leftAsymptote <- c(NA,-10*rightAsymptote[2]) # arbitrary large value
    }
  }
  if((length(na.omit(empDens[1:maxIndex]))<2) &
       (length(na.omit(empDens[maxIndex:length(empDens)]))<2)){
     if(startValues=="BN"|startValues=="SL")
       stop("not enough breaks to estimate asymptotes to log-density")
   }
  if(startValues == "US"){
    svName <- "User Specified"
      if(is.null(ThetaStart)) stop("ThetaStart must be specified")
      if(!is.null(ThetaStart)){
        if(length(ThetaStart)!=4) stop("ThetaStart must contain 4 values")
        if(ThetaStart[2]<=0)
          stop("zeta in ThetaStart must be greater than zero")
        if(ThetaStart[3]<=0)
          stop("delta in ThetaStart must be greater than zero")
        }
    ThetaStart <- c(hyperbPi = ThetaStart[1], zeta = log(ThetaStart[2]),
                      delta = log(ThetaStart[3]), mu = ThetaStart[4])
                   #this gives correct ThetaStart output
                   #when startValues=="US"
  }
  
  if(startValues=="FN"){
    svName <- "Fitted Normal"
    nu <- as.numeric(midpoints[maxIndex])
    mu <- mean(x)
    delta <- sd(x)
    hyperbPi <- (nu - mu)/delta
    zeta <- 1 + hyperbPi^2
    ThetaStart <- c(hyperbPi, log(zeta), log(delta), mu)
  }
  if(startValues=="SL"){
    svName <- "Skew Laplace" 
    llsklp <- function(Theta){
         -sum(log(dskewlap(x, Theta, logPars = TRUE)))
    }
    lSkewAlpha <- log(1/leftAsymptote[2])
    lSkewBeta <- log(abs(1/rightAsymptote[2]))
    skewMu <- midpoints[maxIndex]
    ThetaStart <- c(lSkewAlpha, lSkewBeta, skewMu)
    skewlpOptim <- optim(ThetaStart, llsklp, NULL,
                         method = startMethodSL, hessian = FALSE, ...)
    phi <- 1/exp(skewlpOptim$par[1])
    hyperbGamma <- 1/exp(skewlpOptim$par[2])
    delta <- 0.1 # Take delta to be small
    mu <- skewlpOptim$par[3]
    hyperbPi <- hyperbChangePars(3, 1, c(phi, hyperbGamma, delta, mu))[1]
    zeta <- hyperbChangePars(3, 1, c(phi, hyperbGamma, delta, mu))[2]
    ThetaStart <- c(hyperbPi, log(zeta), log(delta), mu)
  }
  if(startValues=="MoM"){
     svName <- "Method of Moments"
     ThetaStart <- hyperbFitStartMoM(x, startMethodMoM = startMethodMoM, ...)
   }
  if (!(startValues%in%c("US","FN","SL","MoM"))) startValues <- "BN"
  if(startValues=="BN"){
    svName <- "Barndorff-Nielsen 1977"
    phi <- leftAsymptote[2]
    hyperbGamma <- -rightAsymptote[2]
    if(!(is.na(leftAsymptote[1])|is.na(rightAsymptote[1]))){
      mu <- -(leftAsymptote[1]- rightAsymptote[1])/
             (leftAsymptote[2]- rightAsymptote[2])
      intersectionValue <- leftAsymptote[1]+ mu*leftAsymptote[2]
      logModalDens <- log(max(empDens,na.rm=TRUE))
      zeta <- intersectionValue - logModalDens
      if(zeta <= 0) zeta <- 0.1        # This is set arbitrarily
    }else{
      mu <- median(x)
      intersectionValue <- mu
      logModalDens <- log(max(empDens,na.rm=TRUE))
      zeta <- intersectionValue - logModalDens
      if(zeta <= 0) zeta <- 0.1        # This is set arbitrarily
    }
    delta <- zeta/sqrt(phi*hyperbGamma)
    hyperbPi <- hyperbChangePars(3, 1, c(phi, hyperbGamma, delta, mu))[1]      
    ThetaStart <- c(hyperbPi, log(zeta), log(delta), mu)
  }
  names(ThetaStart) <- c("hyperbPi", "lZeta", "lDelta", "mu")
  list(ThetaStart=ThetaStart, breaks = breaks, midpoints = midpoints,
       empDens = empDens, svName=svName)
     
} ## End of hyperbFitStart()

hyperbFitStartMoM <- function(x, startMethodMoM = "Nelder-Mead", ...){
  fun1 <- function(expTheta){
    diff1 <- hyperbMean(expTheta) - mean(x)
  diff1
  }
  fun2 <- function(expTheta){
    diff2 <- hyperbVar(expTheta) - var(x)
  diff2
  }
  fun3 <- function(expTheta){
    diff3 <- hyperbSkew(expTheta) - skewness(x)
    diff3
  }
  fun4 <- function(expTheta){ 
    diff4 <- hyperbKurt(expTheta) - kurtosis(x) 
    diff4
  }
  MoMOptimFun <- function(Theta){
    expTheta <-c(Theta[1],exp(Theta[2]),exp(Theta[3]),Theta[4])
    (fun1(expTheta))^2 + (fun2(expTheta))^2 + 
    (fun3(expTheta))^2 + (fun4(expTheta))^2
  }
  ## Determine starting values for parameters using
  ## Barndorff-Nielsen et al "The Fascination of Sand" in
  ## A Celebration of Statistics pp.78--79
  xi <- sqrt(kurtosis(x)/3)
  chi <- skewness(x)/3  # Ensure 0 <= |chi| < xi < 1
  if(xi >= 1) xi <- 0.999
  adjust <- 0.001
  if(abs(chi) > xi){
    if(xi < 0 ){
    chi <- xi + adjust
    }else{
    chi <- xi - adjust
    }
  }
  hyperbPi <- chi/sqrt(xi^2-chi^2)
  zeta <- 3/xi^2 - 1
  rho <- chi/xi
  delta <- (sqrt(1 + zeta) - 1)* sqrt( 1 - rho^2)
  mu <- mean(x) - delta * hyperbPi * RLambda(zeta, lambda = 1)
  startValuesMoM <- c(hyperbPi, log(zeta), log(delta), mu)
  ## Get Method of Moments estimates
  MoMOptim <- optim(startValuesMoM, MoMOptimFun, method = startMethodMoM, ...)
  ThetaStart <- MoMOptim$par
  ThetaStart
} ## End of hyperbFitStartMoM



