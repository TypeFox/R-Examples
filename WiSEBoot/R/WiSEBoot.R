WiSEBoot <-
function(X, R=100, XParam=NA, TauSq="log", bootDistn="normal", 
                           by.row=FALSE, J0=NA, 
                           wavFam="DaubLeAsymm", wavFil=8, wavBC="periodic"){

  ##Check X is a vector or matrix##
  if(is.matrix(X)!=TRUE && is.atomic(X)!=TRUE){
    stop("X must be a vector or matrix")
  }
  else if((dim(X)[1]==0 || dim(X)[2]==0) && length(X)==0){
    stop("X must have entries")
  }
  else if(mode(X)!="numeric"){
    stop("X must be of type numeric")
  }
  else if(anyNA(X)==TRUE){
    stop("X must not have any missing values.  Please impute first.")
  }
  else if(is.matrix(X)==FALSE && !(log(length(X), base=2)%%1==0) ){
    stop("X must have a length which is a power of 2")
  }
  else if(is.matrix(X)==TRUE && by.row==FALSE && !(log(dim(X)[1], base=2)%%1==0) ){
    stop("The columns of X are the time series.  Time series must have a length which is a power of 2.")
  }
  else if(is.matrix(X)==TRUE && by.row==TRUE && !(log(dim(X)[2], base=2)%%1==0) ){
    stop("The rows of X are the time series.  Time series must have a length which is a power of 2.")
  }


  ##Check by.row is TRUE or FALSE##
  if(!(by.row %in% c(TRUE, FALSE)) ){
    stop("by.row should be logical")
  }


  ####CREATE SOME VARIABLES FOR CONVENIENCE####
  ##If X is a vector, convert to a 1-dimensional matrix, transpose if row contains the data series
  if(is.matrix(X)==FALSE){
    X <- as.matrix(X, ncol=1)
  }else if(by.row==TRUE){ 
    X <- t(X)                
  }
  seriesSamp <- dim(X)[2]     #number of observed data series in the matrix
  J <- log(dim(X)[1], base=2) #power of 2 the length of X is.


  ##Check R is a positive integer##
  if(mode(R)!="numeric" || length(R)>1 || R<=0 || !(R%%1==0) ){
    stop("R must be a positive number.")
  }

  ##Check XParam##
  if(length(XParam)==1 &&!is.na(XParam) ){
    stop("XParam should be missing, a vector, or a matrix.")
  }else if(is.matrix(XParam)!=TRUE && is.atomic(XParam)!=TRUE ){
    stop("XParam should be missing, a vector, or a matrix.")
  }else if(length(XParam)>1 && is.matrix(XParam)==TRUE && dim(XParam)[1]!=2){
    stop("XParam should have 2 rows.")
  }else if(length(XParam)>1 && (anyNA(XParam)==TRUE || mode(XParam)!="numeric")){
    stop("XParam supplied should be numeric with no missing values.")
  }else if(length(XParam)>1 && dim(X)[2]!=1 && (is.matrix(XParam)==FALSE || dim(XParam)[2]!=dim(X)[2])){
    stop("XParam should have the same number of columns as X.")
  }else if(length(XParam)>1 && dim(X)[2]==1 && ( (is.matrix(XParam)==TRUE && dim(XParam)!=c(2,1)) ||
                                                 (is.matrix(XParam)==FALSE && length(XParam)!=2) ) ){ 
    stop("XParam should have 1 column (to match X).")
  }

  ##If XParam is a vector, convert it to a 1-dimensional matrix
  if(is.matrix(XParam)==FALSE){
    XParam <- as.matrix(XParam, ncol=1)
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


  ##Check smoothLevel##
  if(!is.na(J0) && (J0< -1 || !(J0%%1==0))){
    stop("J0 should be NA or an integer which is -1 or greater.")
  }else if(!is.na(J0) && J0>=(J-1)){
    stop("J0 must be less than the power of 2 which is the time series length.")
  }
  if(!is.na(J0)){smoothLevel <- J0+1} else{smoothLevel <- J0}

  ##Check bootVar##
  if(!(TauSq %in% c("log","log10","sqrt","2/5", "1")) ){
    stop("Invalid value for TauSq")
  }

  ####SET SCALE VALUE####
  if(TauSq=="log"){
    vBoot <- log(2^J)
  }else if(TauSq=="log10"){
    vBoot <- log(2^J, base=10)
  }else if(TauSq=="sqrt"){
    vBoot <- sqrt(2^J)
  }else if(TauSq=="1"){
    vBoot <- 1
  }else{
    vBoot <- (2^J)^(2/5)
  }

  ##Check bootDistn##
  if(!( bootDistn %in% c("normal","uniform","exponential","laplace","lognormal","gumbel","t5","t8","t14") )){
    stop("Invalid value for bootDistn.")
  }





  ##Estimate slope (by index), intercept of each data series and keep smoothed residuals, wavelet residuals
  estpIntercept <- matrix(nrow=1, ncol=seriesSamp)             #estimate the population intercept from data
  estpSlope <- matrix(nrow=1, ncol=seriesSamp)                 #estimate the population slope from data

  #keep smoothed residuals and wavelet residuals in arrays
  if(!is.na(smoothLevel)){
    obsDataResidSmooth <- matrix(nrow=2^J, ncol=seriesSamp)
    waveletResiduals <- matrix(nrow=2^J, ncol=seriesSamp)
  }else{
    obsDataResidSmooth <- array(dim=c(2^J, seriesSamp, J)) 
    waveletResiduals <- array(dim=c(2^J, seriesSamp, J)) 
  }
  DataWavelet <- matrix(nrow=2^J, ncol=seriesSamp)                       #keep wavelet coefficients from the input data
  
  #loop for the computation
  for(ss in 1:seriesSamp){
    if(anyNA(XParam)==TRUE){
      obsDataregFit <- lm(X[ , ss]~seq(1, 2^J))  #fit the ith series regression
      estpIntercept[1, ss] <- obsDataregFit$coef[1] #save the intercept for the ith series
      estpSlope[1, ss] <- obsDataregFit$coef[2]     #save the slope for the ith series
      obsDataResid <- as.vector(obsDataregFit$residuals) #keep the ith series residuals
    }else{
      estpIntercept[1, ss] <- XParam[1, ss] #save the intercept for the ith series
      estpSlope[1, ss] <- XParam[2, ss]     #save the slope for the ith series
      obsDataResid <- X[ , ss]              #user supplied a de-trended series
    }
    holdSmoothSeries <- smoothTimeSeries(obsDataResid, wavFam=wavFam, wavFil=wavFil, wavBC=wavBC) #smooth the ss-th series residuals with wavelets
    
    #save the appropriate smoothed residuals
    if(!is.na(smoothLevel)){
      obsDataResidSmooth[ , ss] <- as.matrix(holdSmoothSeries[ , (colnames(holdSmoothSeries) %in% c(paste("J0plusOne",smoothLevel, sep='')))])
      waveletResiduals[ , ss] <- obsDataResid - obsDataResidSmooth[ , ss] #obtain the wavelet residuals

      datawave <- wd(obsDataResidSmooth[ , ss], family=wavFam, filter.number=wavFil, bc=wavBC)
      dw <- accessC(datawave, level=0)
      for( w in 0:(J-1) ){
        dw <- c(dw, accessD(datawave, level=w) )
      }
      DataWavelet[ , ss] <- dw
      remove(dw, datawave)
    }else{
      obsDataResidSmooth[ , ss, ] <- as.matrix(holdSmoothSeries[ , !(colnames(holdSmoothSeries) %in% c("origSeries"))])
      waveletResiduals[ , ss, ] <- obsDataResid - obsDataResidSmooth[ , ss, ] #obtain the wavelet residuals
    }

  }
  remove(holdSmoothSeries, obsDataResid)  #remove extra junk to save memory space




  if(is.na(smoothLevel)){
    bootWaveCoef <- array(dim=c(R, 2^J, seriesSamp, J)) #keep bootstrap wavelet coefficients in an array
    bootIntercept <- array(dim=c(R, seriesSamp, J))     #keep estimated intercepts from bootstrap
    bootSlope <- array(dim=c(R, seriesSamp, J))         #keep estimated slopes from bootstrap
    meanOfMSE <- rep(0, J)                              #evaluation criteria (min) for smoothing level

    for (j in 1:J){
      if(bootDistn=="normal"){
        bootWtMatrix <- sqrt(vBoot)*matrix(rnorm(R*2^J, mean=0, sd=1), nrow=R, ncol=2^J)
      }else if(bootDistn=="uniform"){
        bootWtMatrix <- sqrt(vBoot)*matrix(runif(R*2^J, min=-sqrt(3), max=sqrt(3)), nrow=R, ncol=2^J)
      }else if(bootDistn=="exponential"){
        bootWtMatrix <- sqrt(vBoot)*matrix(rexp(R*2^J)-1, nrow=R, ncol=2^J)
      }else if(bootDistn=="laplace"){
        uniforms <- runif(R*2^J, min=-1/2, max=1/2)
        bootWtMatrix <- sqrt(vBoot)*matrix(-1/sqrt(2)*sign(uniforms)*log(1-2*abs(uniforms)), nrow=R, ncol=2^J)
      }else if(bootDistn=="lognormal"){
        sig <- 1
        mu  <- (log(1/(exp(sig)-1))-sig)/2
        bootWtMatrix <- sqrt(vBoot)*matrix(rlnorm(R*2^J, meanlog=mu, sdlog=sig)-exp(mu+sig/2), nrow=R, ncol=2^J)
      }else if(bootDistn=="gumbel"){
        bootWtMatrix <- sqrt(vBoot)*matrix(rgumbel(R*2^J, scale=sqrt(6/pi^2), location=sqrt(6/pi^2)*digamma(1)), nrow=R, ncol=2^J)
      }else if(bootDistn=="t5"){
        bootWtMatrix <- sqrt(vBoot)*matrix(rt(R*2^J, df=5)*sqrt(3/5), nrow=R, ncol=2^J)
      }else if(bootDistn=="t8"){
        bootWtMatrix <- sqrt(vBoot)*matrix(rt(R*2^J, df=8)*sqrt(6/8), nrow=R, ncol=2^J)
      }else if(bootDistn=="t14"){
        bootWtMatrix <- sqrt(vBoot)*matrix(rt(R*2^J, df=14)*sqrt(12/14), nrow=R, ncol=2^J)
      }
      MSE <- rep(NA, R)
  
      for(r in 1:R){
        bootobsData <- (matrix(1, nrow=2^J, ncol=1) %*% estpIntercept + 
                        matrix(seq(1, 2^J), nrow=2^J, ncol=1) %*% estpSlope + 
                        obsDataResidSmooth[ , , j] + 
                        diag(bootWtMatrix[r, ]) %*% waveletResiduals[ , , j]) 
        bootregFitResidSmooth <- matrix(nrow=2^J, ncol=seriesSamp)

        for(ss in 1:seriesSamp){
          bootregFit <- lm(bootobsData[ , ss]~seq(1, 2^J))   #estimate the slope and intercept in the i-th boot data
          bootIntercept[r, ss, j] <-  bootregFit$coef[1]     #save the estimated bootstrap intercept
          bootSlope[r, ss, j] <- bootregFit$coef[2]          #save the estimated bootstrap slope
          bootregFitResid <- as.vector(bootregFit$residuals) #use de-trended series in wavelet smoothing for bootstrap
    
          bootregFitResidWaveDecomp <- wd(bootregFitResid, filter.number=wavFil, family=wavFam, type="wavelet", bc=wavBC)
          for(k in (J-1):(J-j)){
            bootregFitResidWaveDecomp <- putD(bootregFitResidWaveDecomp, level=k, rep(0, 2^k))  #put zeroes into the wavelet object, at the correct level
          }
  
          # wave wavelet coefficients from the smoothed object
          bwv <- accessC(bootregFitResidWaveDecomp, level=0)
          for(k in 0:(J-1)){
            bwv <- c(bwv, accessD(bootregFitResidWaveDecomp, level=k))
          }
          bootWaveCoef[r, , ss, j] <- bwv

          bootregFitResidSmooth[ , ss] <- wr(bootregFitResidWaveDecomp) #construct the smoothed bootstrap sample
        }
        diffMat <- (X - (matrix(1, nrow=2^J, ncol=1) %*% bootIntercept[r, , j] + 
                         matrix(seq(1, 2^J), nrow=2^J, ncol=1) %*% bootSlope[r, , j] + 
                            bootregFitResidSmooth) )
        MSE[r] <- sum(diag( t(diffMat) %*% diffMat ))/(2^J*seriesSamp)
      }

      meanOfMSE[j] <- sum(MSE)/length(MSE)
    }

    MSECriteria <- as.matrix(cbind(J0plusOne=seq(J-1,0), meanOfMSE))
    bootIntercept <- bootIntercept[ , , which.min(MSECriteria[ , 2])]
    bootSlope <- bootSlope[ , , which.min(MSECriteria[ , 2])]
    bootWaveCoef <- bootWaveCoef[ , , , which.min(MSECriteria[ , 2])]

    for (ss in 1:seriesSamp){
      datawave <- wd(obsDataResidSmooth[ , ss, which.min(MSECriteria[ , 2])], family=wavFam, 
                       filter.number=wavFil, bc=wavBC)
      dw <- accessC(datawave, level=0)
      for( w in 0:(J-1) ){
        dw <- c(dw, accessD(datawave, level=w) )
      }
      DataWavelet[ , ss] <- dw
    }


  }else{
    bootIntercept <- matrix(nrow=R, ncol=seriesSamp)
    bootSlope <- matrix(nrow=R, ncol=seriesSamp)
    bootWaveCoef <- array(dim=c(R, 2^J, seriesSamp) )
    if(bootDistn=="normal"){
      bootWtMatrix <- sqrt(vBoot)*matrix(rnorm(R*2^J, mean=0, sd=1), nrow=R, ncol=2^J)
    }else if(bootDistn=="uniform"){
      bootWtMatrix <- sqrt(vBoot)*matrix(runif(R*2^J, min=-sqrt(3), max=sqrt(3)), nrow=R, ncol=2^J)
    }else if(bootDistn=="exponential"){
      bootWtMatrix <- sqrt(vBoot)*matrix(rexp(R*2^J)-1, nrow=R, ncol=2^J)
    }else if(bootDistn=="laplace"){
      uniforms <- runif(R*2^J, min=-1/2, max=1/2)
      bootWtMatrix <- sqrt(vBoot)*matrix(-1/sqrt(2)*sign(uniforms)*log(1-2*abs(uniforms)), nrow=R, ncol=2^J)
    }else if(bootDistn=="lognormal"){
      bootWtMatrix <- sqrt(vBoot)*matrix(rlnorm(R*2^J, meanlog=1-sig^2/2, sdlog=sig)-exp(1), nrow=R, ncol=2^J)
    }else if(bootDistn=="gumbel"){
      bootWtMatrix <- sqrt(vBoot)*matrix(rgumbel(R*2^J, scale=sqrt(6/pi^2), location=sqrt(6/pi^2)*digamma(1)), nrow=R, ncol=2^J)
    }else if(bootDistn=="t5"){
      bootWtMatrix <- sqrt(vBoot)*matrix(rt(R*2^J, df=5)*sqrt(3/5), nrow=R, ncol=2^J)
    }else if(bootDistn=="t8"){
      bootWtMatrix <- sqrt(vBoot)*matrix(rt(R*2^J, df=8)*sqrt(6/8), nrow=R, ncol=2^J)
    }else if(bootDistn=="t14"){
      bootWtMatrix <- sqrt(vBoot)*matrix(rt(R*2^J, df=14)*sqrt(12/14), nrow=R, ncol=2^J)
    }
    MSE <- rep(0, R)

    for(r in 1:R){
      bootobsData <- (matrix(1, nrow=2^J, ncol=1) %*% estpIntercept + 
                        matrix(seq(1, 2^J), nrow=2^J, ncol=1) %*% estpSlope + 
                        obsDataResidSmooth + 
                        diag(bootWtMatrix[r, ]) %*% waveletResiduals) 
        bootregFitResidSmooth <- matrix(nrow=2^J, ncol=seriesSamp)

        for(ss in 1:seriesSamp){
          bootregFit <- lm(bootobsData[ , ss]~seq(1, 2^J))   #estimate the slope and intercept in the i-th boot data
          bootIntercept[r, ss] <-  bootregFit$coef[1]        #save the estimated bootstrap intercept
          bootSlope[r, ss] <- bootregFit$coef[2]             #save the estimated bootstrap slope
          bootregFitResid <- as.vector(bootregFit$residuals) #use de-trended series in wavelet smoothing for bootstrap
    
          bootregFitResidWaveDecomp <- wd(bootregFitResid, filter.number=wavFil, family=wavFam, type="wavelet", bc=wavBC)
          for(k in (J-1):smoothLevel){
            bootregFitResidWaveDecomp <- putD(bootregFitResidWaveDecomp, level=k, rep(0, 2^k))  #put zeroes into the wavelet object, at the correct level
          }

          # wave wavelet coefficients from the smoothed object
          bwv <- accessC(bootregFitResidWaveDecomp, level=0)
          for(k in 0:(J-1)){
            bwv <- c(bwv, accessD(bootregFitResidWaveDecomp, level=k))
          }
          bootWaveCoef[r, , ss] <- bwv

        bootregFitResidSmooth[ , ss] <- wr(bootregFitResidWaveDecomp) #construct the smoothed bootstrap sample
        }

        diffMat <- (X - (matrix(1, nrow=2^J, ncol=1) %*% bootIntercept[r, ] + 
                         matrix(seq(1, 2^J), nrow=2^J, ncol=1) %*% bootSlope[r, ] + 
                            bootregFitResidSmooth) )
        MSE[r] <- sum(diag( t(diffMat) %*% diffMat ))/(2^J*seriesSamp)
    }
    meanOfMSE <- sum(MSE)/length(MSE)
    MSECriteria <- as.matrix(cbind(J0plusOne=smoothLevel, meanOfMSE))
  }


  structure(invisible(list( MSECriteria=MSECriteria, BootIntercept=bootIntercept, 
                 BootSlope=bootSlope, BootWavelet=bootWaveCoef, DataWavelet=DataWavelet,
                 wavFam=wavFam, wavFil=wavFil, wavBC=wavBC, TauSq=TauSq, bootDistn=bootDistn, by.row=by.row,
                 XParam=t(matrix(cbind(estpIntercept, estpSlope), ncol=2)) )), 
            class="WiSEBoot")

}
