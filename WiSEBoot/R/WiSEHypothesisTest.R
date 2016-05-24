WiSEHypothesisTest <-
function(X, Y, J0, R=100, popParam=c(0,1), XParam=c(NA,NA), YParam=c(NA, NA),
                               TauSq="log", bootDistn="normal", 
                               wavFam="DaubLeAsymm", wavFil=8, wavBC="periodic", plot=TRUE, ...){

  ##Check Y and X##
  if(is.atomic(X)!=TRUE || is.atomic(Y)!=TRUE){
    stop("Y and X must be a vectors.")
  }
  else if( (is.matrix(Y)==TRUE && dim(Y)[2]>1) || (is.matrix(X)==TRUE && dim(X)[2]>1) ){
    stop("Y and X supplied as matrices must each have only one column.")
  }
  else if(length(X)==0 || length(Y)==0){
    stop("Y and X must have entries.")
  }
  else if(mode(X)!="numeric" || mode(Y)!="numeric"){
    stop("Y and X must be of type numeric.")
  }
  else if(anyNA(X)==TRUE || anyNA(Y)==TRUE){
    stop("Y and X must not have any missing values.")
  }
  else if(length(X)!=length(Y)){
    stop("Y and X must have the same number of observations.")
  }
  else if(!(log(length(X), base=2)%%1==0) ){
    stop("Y and X must have a length which is a power of 2")
  }
  else if(length(X)<=2){
    stop("This method is not useful for vectors of length less than 4.")
  }

  ##If Y or X is a vector, convert to a 1-dimensional matrix##
  if(is.matrix(X)==FALSE){ X <- as.matrix(X, ncol=1) }
  if(is.matrix(Y)==FALSE){ Y <- as.matrix(Y, ncol=1) }

  J <- log(dim(X)[1], base=2) #power of 2 the length of data is.


  ##Check smooth level (J0) is between 0 and J-2##
  if( J0 < 0 || J0 > J-2 || !(J0%%1==0) || length(J0)>1 ){
    stop("The threshold level, J0, should be an integer between 0 and J-2.")
  }else if(is.na(J0)){
    stop("The threshold level must be set.  To select a threshold level use the 
          WiSE.bootstrap function.")
  }
  

  ##Check R is a positive integer##
  if(mode(R)!="numeric" || length(R)>1 || R<=0 || !(R%%1==0) ){
    stop("R must be a positive number.")
  }

  ##Check boundary condition##
  if(!(wavBC %in% c("periodic","symmetric")) ){
    stop("wavBC may only take values periodic or symmetric.")
  }

  ##Check wavelet filter and family##
  if(!( (wavFam=="DaubLeAsymm" && wavFil %in% seq(4,10)) || 
        (wavFam=="DaubExPhase" && wavFil %in% seq(1,10)) )){
    stop("wavFam and wavFil combination not allowed in wavethresh.")
  }


  ##Check bootVar##
  if(!(TauSq %in% c("log","log10","sqrt","2/5")) ){
    stop("Invalid value for TauSq")
  }

  ####SET SCALE VALUE####
  if(TauSq=="log"){
    vBoot <- log(2^J)
  }else if(TauSq=="log10"){
    vBoot <- log(2^J, base=10)
  }else if(TauSq=="sqrt"){
    vBoot <- sqrt(2^J)
  }else if (TauSq=="1"){
    vBoot <- 1
  }else{
    vBoot <- (2^J)^(2/5)
  }


  ##Check bootDistn##
  if(!( bootDistn %in% c("normal","uniform","exponential","laplace","lognormal","gumbel","t5","t8","t14") )){
    stop("Invalid value for bootDistn.")
  }

  ##Check popParam##
  if(is.atomic(popParam)!=TRUE || is.null(dim(popParam))!=TRUE){
    stop("popParam should be a vector.")
  }
  else if(anyNA(popParam) || length(popParam)!=2 || mode(popParam)!="numeric"){
    stop("popParam should contain 2 numbers.")
  }

  ##Check XParam and YParam##
  if(is.atomic(XParam)!=TRUE || is.null(dim(XParam))!=TRUE){
    stop("XParam should be a vector.")
  }
  else if(length(XParam)!=2){
    stop("XParam should be a vector of length 2.")
  }
  else if(sum(is.na(XParam))==1){
    stop("XParam should contain 2 numbers of be completely missing.")
  }
  else if(sum(is.na(XParam))==0 && mode(XParam)!="numeric" ){
    stop("XParam should contain 2 numbers or be completely missing.")
  }
  if(is.atomic(YParam)!=TRUE || is.null(dim(YParam))!=TRUE){
    stop("YParam should be a vector.")
  }
  else if(length(XParam)!=2){
    stop("XParam should be a vector of length 2.")
  }
  else if(sum(is.na(YParam))==1){
    stop("YParam should contain 2 numbers of be completely missing.")
  }
  else if(sum(is.na(YParam))==0 && mode(YParam)!="numeric" ){
    stop("YParam should contain 2 numbers or be completely missing.")
  }


  ##Check plot##
  if(!(plot %in% c(TRUE, FALSE))){
    stop("plot is logical.")
  }
  dots <- list(...)




  ##Estimate least squares intercept and slope (by index) for each series
  if(anyNA(YParam)==TRUE){
    YregFit <- lm(Y~seq(1, 2^J))  
    estYIntercept <- YregFit$coef[1]
    estYSlope <- YregFit$coef[2]
    YResid <- as.vector(YregFit$residuals)
  }else{
    estYIntercept <- YParam[1]
    estYSlope <- YParam[2]
    YResid <- Y
  }

  if(anyNA(XParam)==TRUE){
    XregFit <- lm(X~seq(1, 2^J))
    estXIntercept <- XregFit$coef[1]
    estXSlope <- XregFit$coef[2]
    XResid <- as.vector(XregFit$residuals)
  }else{
    estXIntercept <- XParam[1]
    estXSlope <- XParam[2]
    XResid <- X
  }



  ##Data wavelet decomposition
  YWave <- wd(YResid, family=wavFam, filter.number=wavFil, bc=wavBC)
  XWave <- wd(XResid, family=wavFam, filter.number=wavFil, bc=wavBC)
  YScalingCoef <- accessC(YWave, level=0)      #scaling coefficient from the Y (to be kept in the wavelet)

  

  ##Original wavelet coefficient estimates from the data
  estYWavelet <- rep(NA, 2^(J0+1)-1)   #keep the non-thresholded mother wavelet coefficients for Y
  estXWavelet <- rep(NA, 2^(J0+1)-1)   #keep the non-thresholded mother wavelet coefficients for X
  for(j in 0:J0){
    estYWavelet[(2^j):(2^(j+1)-1)] <- accessD(YWave, level=j)
    estXWavelet[(2^j):(2^(j+1)-1)] <- accessD(XWave, level=j)
  }
  meanYWavelet <- sum(estYWavelet)/length(estYWavelet)  #mean of the coarse coefficients for Y data
  meanXWavelet <- sum(estXWavelet)/length(estXWavelet)  #mean of the coarse coefficients for X data



  ##Parameter estimates from the data
  dataSlope <- sum((estYWavelet - meanYWavelet)*(estXWavelet - meanXWavelet))/sum((estXWavelet - meanXWavelet)^2)
  dataIntercept <- meanYWavelet - dataSlope*meanXWavelet



  ##Smooth the X series to the selected level
  holdXSmoothWave <- XWave
  holdYSmoothWave <- YWave
  for(j in (J0+1):(J-1)){
    holdXSmoothWave <- putD(holdXSmoothWave, level=j, rep(0, 2^j))
    holdYSmoothWave <- putD(holdYSmoothWave, level=j, rep(0, 2^j))
  }
  XResidSmooth <- wr(holdXSmoothWave)
  XwaveletResiduals <- XResid - XResidSmooth
  YResidSmooth <- wr(holdYSmoothWave)
  YwaveletResiduals <- YResid - YResidSmooth






  ##Save bootstrap estimated filter wavelet coefficients from each series
  bootEstYWavelet <- matrix(nrow=R, ncol=length(estYWavelet))
  bootEstXWavelet <- matrix(nrow=R, ncol=length(estYWavelet))
  bootSlope <- rep(NA, R)
  bootIntercept <- rep(NA, R)

  ##Set the bootstrap weight matrix based upon specified distribution
  if(bootDistn=="normal"){
    bootWtMatrix <- sqrt(vBoot)*matrix(rnorm(R*2^J, mean=0, sd=1), nrow=R, ncol=2^J)
    bootWtMatrix2 <- sqrt(vBoot)*matrix(rnorm(R*2^J, mean=0, sd=1), nrow=R, ncol=2^J)
  }else if(bootDistn=="uniform"){
    bootWtMatrix <- sqrt(vBoot)*matrix(runif(R*2^J, min=-sqrt(3), max=sqrt(3)), nrow=R, ncol=2^J)
    bootWtMatrix2 <- sqrt(vBoot)*matrix(runif(R*2^J, min=-sqrt(3), max=sqrt(3)), nrow=R, ncol=2^J)
  }else if(bootDistn=="exponential"){
    bootWtMatrix <- sqrt(vBoot)*matrix(rexp(R*2^J)-1, nrow=R, ncol=2^J)
    bootWtMatrix2 <- sqrt(vBoot)*matrix(rexp(R*2^J)-1, nrow=R, ncol=2^J)
  }else if(bootDistn=="laplace"){
    uniforms <- runif(R*2^J, min=-1/2, max=1/2)
    bootWtMatrix <- sqrt(vBoot)*matrix(-1/sqrt(2)*sign(uniforms)*log(1-2*abs(uniforms)), nrow=R, ncol=2^J)
    uniforms2 <- runif(R*2^J, min=-1/2, max=1/2)
    bootWtMatrix2 <- sqrt(vBoot)*matrix(-1/sqrt(2)*sign(uniforms)*log(1-2*abs(uniforms)), nrow=R, ncol=2^J)
  }else if(bootDistn=="lognormal"){
    sig <- 1
    mu  <- (log(1/(exp(sig)-1))-sig)/2
    bootWtMatrix <- sqrt(vBoot)*matrix(rlnorm(R*2^J, meanlog=mu, sdlog=sig)-exp(mu+sig/2), nrow=R, ncol=2^J)
    bootWtMatrix2 <- sqrt(vBoot)*matrix(rlnorm(R*2^J, meanlog=mu, sdlog=sig)-exp(mu+sig/2), nrow=R, ncol=2^J)
  }else if(bootDistn=="gumbel"){
    bootWtMatrix <- sqrt(vBoot)*matrix(rgumbel(R*2^J, scale=sqrt(6/pi^2), location=sqrt(6/pi^2)*digamma(1)), nrow=R, ncol=2^J)
    bootWtMatrix2 <- sqrt(vBoot)*matrix(rgumbel(R*2^J, scale=sqrt(6/pi^2), location=sqrt(6/pi^2)*digamma(1)), nrow=R, ncol=2^J)
  }else if(bootDistn=="t5"){
    bootWtMatrix <- sqrt(vBoot)*matrix(rt(R*2^J, df=5)*sqrt(3/5), nrow=R, ncol=2^J)
    bootWtMatrix2 <- sqrt(vBoot)*matrix(rt(R*2^J, df=5)*sqrt(3/5), nrow=R, ncol=2^J)
  }else if(bootDistn=="t8"){
    bootWtMatrix <- sqrt(vBoot)*matrix(rt(R*2^J, df=8)*sqrt(6/8), nrow=R, ncol=2^J)
    bootWtMatrix2 <- sqrt(vBoot)*matrix(rt(R*2^J, df=8)*sqrt(6/8), nrow=R, ncol=2^J)
  }else if(bootDistn=="t14"){
    bootWtMatrix <- sqrt(vBoot)*matrix(rt(R*2^J, df=14)*sqrt(12/14), nrow=R, ncol=2^J)
    bootWtMatrix2 <- sqrt(vBoot)*matrix(rt(R*2^J, df=14)*sqrt(12/14), nrow=R, ncol=2^J)
  }

  for(r in 1:R){
    bootXData <- (matrix(1, nrow=2^J, ncol=1) %*% estXIntercept + 
                         matrix(seq(1, 2^J), nrow=2^J, ncol=1) %*% estXSlope + 
                         XResidSmooth + diag(bootWtMatrix[r, ]) %*% XwaveletResiduals) 

    bootXregFit <- lm(bootXData~seq(1, 2^J))
    bootXregFitResid <- as.vector(bootXregFit$residuals)    
    bootXregFitResidWave <- wd(bootXregFitResid, filter.number=wavFil, family=wavFam, type="wavelet", bc=wavBC)

    #Smooth the X residuals
    for(j in (J0+1):(J-1)){
      bootXregFitResidWave <- putD(bootXregFitResidWave, level=j, rep(0, 2^j))
    }
    
    #get the smoothed series so we can bootstrap Y
    holdYSmoothWave <- putC(bootXregFitResidWave, level=0, YScalingCoef)
    for(j in 0:J0){
      currentXCoef <- accessD(holdYSmoothWave, level=j)
      holdYSmoothWave <- putD(holdYSmoothWave, level=j, (popParam[1] + popParam[2]*currentXCoef))
    }
    YResidSmooth <- wr(holdYSmoothWave)
    


    bootYData <- (matrix(1, nrow=2^J, ncol=1) %*% estYIntercept + 
                         matrix(seq(1, 2^J), nrow=2^J, ncol=1) %*% estYSlope + 
                         YResidSmooth + diag(bootWtMatrix2[r, ]) %*% YwaveletResiduals) 
    bootYregFit <- lm(bootYData~seq(1, 2^J))
    bootYregFitResid <- as.vector(bootYregFit$residuals) 
    bootYregFitResidWave <- wd(bootYregFitResid, filter.number=wavFil, family=wavFam, type="wavelet", bc=wavBC)

    for(j in 0:J0){
      bootEstYWavelet[r, (2^j):(2^(j+1)-1)] <- accessD(bootYregFitResidWave, level=j)
      bootEstXWavelet[r, (2^j):(2^(j+1)-1)] <- accessD(bootXregFitResidWave, level=j)
    }


    ##Calculate bootstrap parameter estimates
    meanBootYWavelet <- sum(bootEstYWavelet[r, ])/length(bootEstYWavelet[r, ])
    meanBootXWavelet <- sum(bootEstXWavelet[r, ])/length(bootEstXWavelet[r, ])
    bootSlope[r] <- (sum((bootEstYWavelet[r, ]-meanBootYWavelet)*(bootEstXWavelet[r, ]-meanBootXWavelet))/
                     sum((bootEstXWavelet[r, ]-meanBootXWavelet)^2) )
    bootIntercept[r] <- meanBootYWavelet - bootSlope[r]*meanBootXWavelet
  }


  ##Calculate a p-value associated with asymptotic normality
  invS <- solve(matrix(cbind(c(var(bootIntercept), cov(bootIntercept, bootSlope)),
                             c(cov(bootIntercept, bootSlope), var(bootSlope))), ncol=2))
  diffMat <- matrix(c(dataIntercept-popParam[1], dataSlope-popParam[2]), ncol=1)
  Tsq <- vBoot * t(diffMat) %*% invS %*% diffMat
  numbCoef <- length(estYWavelet)  #the number of coefficients used in the linear regression to obtain alpha, beta

  BootParam <- matrix(cbind(bootIntercept, bootSlope), ncol=2)
  BootParam <- t(BootParam)
  bootDiffMat <- BootParam - popParam  
  bootTsq <- vBoot*t(bootDiffMat) %*% invS %*% bootDiffMat
  bootTsq <- as.vector(diag(bootTsq))


  APValue <- pf((numbCoef-2)/2/(numbCoef-1)*Tsq, df1=2, df2=numbCoef-2, lower.tail=FALSE) 
  BPValue <- 1 - (rank(c(Tsq, bootTsq), ties.method="average")[1] - 1)/R 

  if(plot==TRUE){
    dots$x <- bootIntercept
    dots$y <- bootSlope
    if(is.null(dots$main)==TRUE){dots$main <- paste("Bootstrap Parameter Estimates (under Null), p-value=", prettyNum(BPValue, format="g", digits=3))}
    if(is.null(dots$xlab)==TRUE){dots$xlab <- expression(paste("Intercept (", alpha, ")", sep=''))}
    if(is.null(dots$ylab)==TRUE){dots$ylab <- expression(paste("Slope (", beta, ")", sep=''))}
    if(is.null(dots$xlim)==TRUE){dots$xlim <- c(min(bootIntercept, dataIntercept), max(bootIntercept, dataIntercept))}
    if(is.null(dots$ylim)==TRUE){dots$ylim <- c(min(bootSlope, dataSlope), max(bootSlope, dataSlope))}
    if(is.null(dots$pch)==TRUE){dots$pch <- 20}
    do.call("plot", dots)
    abline(v=popParam[1], col="lightgray")
    abline(h=popParam[2], col="lightgray")
    points(dataIntercept, dataSlope, pch=20, cex=2, col="red")
  }

  structure(invisible(list(AsymptoticPValue=APValue, BootstrapPValue=BPValue,
                           dataSlope=dataSlope, dataIntercept=dataIntercept, 
                           bootSlope=bootSlope, bootIntercept=bootIntercept,
                           YWavelet=estYWavelet, XWavelet=estXWavelet,
                           bootYWavelet=bootEstYWavelet, bootXWavelet=bootEstXWavelet)), class="WiSEHypothesisTest")


}