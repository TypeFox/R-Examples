# Reference:
# Kourentzes N., Petropoulos F. and Trapero J.R. (2014) 
# Improving forecasting by estimating time series structural 
# components across multiple frequencies. 
# International Journal of Forecasting, 30 (2), 291-302
#
# Nikolaos Kourentzes & Fotios Petropoulos (2014)

library(forecast, quietly=TRUE)   # ETS method
library(parallel, quietly=TRUE)   # Needed for parallel
# library(miscTools, quietly=TRUE)  # Needed for medians by columns - Replaced colMedians(X) by apply(X,2,"median")

#-------------------------------------------------
mapa <- function(y, ppy=NULL, fh=ppy, ifh=1, minimumAL=1, maximumAL=ppy, 
	comb="mean", paral=0, display=0, outplot=1, hybrid=TRUE, model="ZZZ", 
  conf.lvl=NULL)
{
# Wrapper to estimate and produce MAPA in- and out-of-sample forecasts
# Uses mapaest and mapafor
#  
# Inputs:
#   y           = In sample observations of a time series (vector)
#                 If y == "paper" then it prints paper reference
#   ppy         = Periods in a season of the time series at the sampled frequency.
#                 If y is a ts object then this is taken from its frequency,
#                 unless overriden. 
#   fh          = Forecast horizon. Default = ppy
#   ifh         = In-sample forecast horizon. Default = 1
#   minimumAL   = Lowest aggregation level to use. Default = 1
#   maximumAL   = Highest aggregation level to use. Default = ppy, maximumAL>1
#   comb        = Combination operator. One of "mean" or "median". Default is "mean"
#   paral       = Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 
#                 2 = yes and initialise cluster. Default is 0.
#   display     = Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
#   outplot     = Provide output plot. 0 = no; 1 = yes. Default is 1.
#   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE
#                 If minimumAL > 1 then the minimumAL ETS forecasts are used.  
#   model       = Allow only that type of ETS at each aggregation level. This follows similar
#                 coding to the ets function. The first letter refers to the error 
#                 type ("A", "M" or "Z"); the second letter refers to the trend 
#                 type ("N","A","Ad","M","Md" or "Z"); and the third letter 
#                 refers to the season type ("N","A","M" or "Z"). The letters mean: 
#                 "N"=none, "A"=additive, "M"=multiplicative and "Z"=automatically selected. 
#                 A "d" for trend implies damped. By default model="ZZZ". If due to sample
#                 limitation ETS cannot be calculated at an aggregation level for the selected
#                 model, then no estimation is done for that specific level. For aggregation 
#                 levels that seasonality becomes 1 then a non-seasonal model is estimated.
#   conf.lvl    = Vector of confidence level for prediction intervals. Values must be (0,1). 
#                 If conf.lvl == NULL then no intervals are calculated. For example to get 
#                 the intervals for 80% and 95% use conf.lvl=c(0.8,0.95). 
#
# Output:
#   out$infor   = In-sample forecasts
#   out$outfor  = Out-of-sample forecasts
#   out$PI      = Prediction intervals for given confidence levels
#   out$MSE     = In-sample MSE error
#   out$MAE     = In-sample MAE error
  
  # Paper info
  if (!is.numeric(y)){
    writeLines("Paper reference: ")
    writeLines("Kourentzes N., Petropoulos F. and Trapero J.R. (2014)")
    writeLines("Improving forecasting by estimating time series structural components")
    writeLines(paste("across multiple frequencies. International Journal of Forecasting,", 
                     " 30 (2), 291-302.",sep=""))
    return(invisible())
  }
  
  # Get ppy, fh and maximumAL
  if (is.null(ppy)){
    if (class(y)=="ts"){
      ppy <- frequency(y)
      if (is.null(fh)){fh <- ppy}
	  if (is.null(maximumAL)){maximumAL <- ppy}
    } else {
      stop(paste("Input ppy is not given and y input not ts class.",
                 "Please provide the periods in a season of the time series",
                 "at the sampled frequency."))
    }
  }
  
  # Estimate MAPA
  mapafit <- mapaest(y=y, ppy=ppy, minimumAL=minimumAL, 
                     maximumAL=maximumAL, paral=paral, display=display, 
                     outplot=outplot, model=model)
  
  # Produce in- and out-of-sample forecasts
  out <- mapafor(y=y, mapafit=mapafit, fh=fh, ifh=ifh, comb=comb, 
                 outplot=outplot, hybrid=hybrid, conf.lvl=conf.lvl)
  
  return(out)
  
}

#-------------------------------------------------
mapasimple <- function(y, ppy=NULL, fh=ppy, minimumAL=1, maximumAL=ppy, comb="mean", 
	output="forecast", paral=0, display=0, outplot=1, hybrid=TRUE, model="ZZZ") 
{
# MAPA estimation and forecast
#  
# Inputs:
#   y           = In sample observations of a time series (vector)
#                 If y == "paper" then it prints paper reference
#   ppy         = Periods in a season of the time series at the sampled frequency.
#                 If y is a ts object then this is taken from its frequency,
#                 unless overriden. 
#   fh          = Forecast horizon. Default = ppy
#   minimumAL   = Lowest aggregation level to use. Default = 1, maximumAL>1
#   maximumAL   = Highest aggregation level to use. Default = ppy
#   comb        = Combination operator. One of "mean" or "median". Default is "mean"
#   output      = Type of output. One of "forecast" or "all". Default is "forecast"
#                 If output="all", both forecasts and components estimates per aggregation
#                 level are provided.
#   paral       = Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 
#                 2 = yes and initialise cluster. Default is 0.
#   display     = Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
#   outplot     = Provide output plot. 0 = no; 1 = time series and forecast only;
#                 2 = time series, forecasts and components. For the components the rainbow 
#                 colouring scheme is used. Red is aggregation level 1, followed by yellow, 
#                 green, cyan, blue and magenta for the higher aggregation levels. Default is 1. 
#   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE  
#                 If minimumAL > 1 then the minimumAL ETS forecasts are used.  
#   model       = Allow only that type of ETS at each aggregation level. This follows similar
#                 coding to the ets function. The first letter refers to the error 
#                 type ("A", "M" or "Z"); the second letter refers to the trend 
#                 type ("N","A","Ad","M","Md" or "Z"); and the third letter 
#                 refers to the season type ("N","A","M" or "Z"). The letters mean: 
#                 "N"=none, "A"=additive, "M"=multiplicative and "Z"=automatically selected. 
#                 A "d" for trend implies damped. By default model="ZZZ". If due to sample
#                 limitation ETS cannot be calculated at an aggregation level for the selected
#                 model, then no estimation is done for that specific level. For aggregation 
#                 levels that seasonality becomes 1 then a non-seasonal model is estimated.
#
# Output:
#   forecasts   = Vector with forecasts
#   components  = array with MAPA components

  # Paper info
  if (!is.numeric(y)){
    writeLines("Paper reference: ")
    writeLines("Kourentzes N., Petropoulos F. and Trapero J.R. (2014)")
    writeLines("Improving forecasting by estimating time series structural components")
    writeLines(paste("across multiple frequencies. International Journal of Forecasting,", 
		                 " 30 (2), 291-302.",sep=""))
    return(invisible())
  }  
  
  # Get ppy, fh and maximumAL
  if (is.null(ppy)){
    if (class(y)=="ts"){
      ppy <- frequency(y)
      if (is.null(fh)){fh <- ppy}
	  if (is.null(maximumAL)){maximumAL <- ppy}
    } else {
    stop(paste("Input ppy is not given and y input not ts class.",
               "Please provide the periods in a season of the time series",
               "at the sampled frequency."))
    }
  }
  
  # Make sure that maximumAL > 1
  if (maximumAL == 1){
    maximumAL <- maximumAL + 1
  }
  
  if (minimumAL>=maximumAL){
    stop("maximumAL must be larger than minimumAL")
  }
  
  # Setup parallel processing if required
  if (paral == 2){
    crs <- detectCores()
    cl <- makeCluster(getOption("cl.cores", crs))
    writeLines(paste("Running with", crs, 'cores'))
  }
  
  observations <- length(y) # number of observations for the in-sample data
  FCs <- array(0, c((maximumAL-minimumAL+1), 4, fh)) # the forecasts and the forecasts of the components 
									  # will be saved here
   
  # Aggregation and estimation 
  if (paral != 0){  # Parallel run
    FCs_par <- clusterApplyLB(cl, 1:(maximumAL-minimumAL+1), mapasimple.loop, 
      y=y, minimumAL=minimumAL, maximumAL=maximumAL, 
      observations=observations, ppy=ppy, display=display, fh=fh, model=model)  
  } else {          # Serial run
    FCs_par <- vector("list", (maximumAL-minimumAL+1))
    for (i in minimumAL:maximumAL){
      FCs_par[[i]] <- mapasimple.loop(i, y, minimumAL, maximumAL, observations, 
        ppy, display, fh, model)
    }
  }
  
  if (paral == 2){
    # Stop parallel processing
    stopCluster(cl)
  }
  
  # Reshape parallel output
  FCs_par <- do.call(rbind, FCs_par)
  
  FCs <- array(0, c((maximumAL-minimumAL+1), 4, fh),dimnames=list(paste("AL",minimumAL:maximumAL,sep=""),
    c("ETS","Level","Trend","Season"),paste("t+",1:fh,sep=""))) # the forecasts and the forecasts 
															                              # of the components will be saved here
  for (f in 1:fh){  
    FCs[, , f] <- t(array(FCs_par[,f],c(4,(maximumAL-minimumAL+1))))
  }
  
  # Check whether all aggregation levels were calculated for given model
  # due to sample size
  names.r <- rownames(FCs)
  names.c <- colnames(FCs)
  AL.idx <- is.na(FCs[,1,1])
  maximumAL <- maximumAL - sum(AL.idx)
  FCs <- FCs[!AL.idx, , ]
  names.r <- names.r[!AL.idx]
  # Check if dimensionality of FCs remains
  if (length(dim(FCs))==2){
    # With fh = 1 the line above removes the third dimension. Add it again
    # Another case is that there is only a single aggregation level left
    FCs <- array(FCs,c(maximumAL,4,fh),dimnames=list(names.r,names.c,paste("t+",1:fh,sep="")))
  }
  
  # MAPA combination
  combres <- mapacomb(minimumAL,maximumAL,ppy,FCs,comb,observations)
  forecasts <- combres[[1]]
  perm_levels <- combres[[2]]
  perm_seas <- combres[[3]]
  
  # Calculate hybrid model
  if (hybrid==TRUE){
    forecasts <- (FCs[1,1,] + forecasts)/2
  }
  
  # Plot output
  mapaplot(outplot,FCs,maximumAL,perm_levels,perm_seas,observations,y,forecasts,fh,comb)
  
  # Construct output
  if (output=="forecast"){
    return(forecasts)
  } else {
    return(list(forecast=forecasts,components=FCs))
  }
  
}

#-------------------------------------------------
mapafor <- function(y, mapafit, fh=-1, ifh=1, comb="mean", outplot=1, hybrid=TRUE,
                    conf.lvl=NULL) {
# MAPA in- and out-of-sample forecast
# 
# Inputs:
#   y           = In sample observations of a time series (vector)
#   mapafit     = Fitted MAPA model (from mapaest)
#   fh          = Forecast horizon. Default = ppy
#   ifh         = In-sample forecast horizon. Default = 1
#   comb        = Combination operator. One of "mean" or "median". Default is "mean"
#   outplot     = Provide output plot. 0 = no; 1 = yes. Default is 1. 
#   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE
#                 If minimumAL > 1 then the minimumAL ETS forecasts are used.  
#   conf.lvl    = Vector of confidence level for prediction intervals. Values must be (0,1). 
#                 If conf.lvl == NULL then no intervals are calculated. For example to get 
#                 the intervals for 80% and 95% use conf.lvl=c(0.8,0.95). 
#
# Output:
#   out$infor   = In-sample forecasts.
#   out$outfor  = Out-of-sample forecasts.
#   out$PI      = Prediction intervals for given confidence levels.
#   out$MSE     = In-sample MSE error.
#   out$MAE     = In-sample MAE error.
  
  observations <- length(y) # number of observations for the in-sample data
  
  # Get settings from mapafit
  ALs <- as.numeric(mapafit[mapafit[,19]==TRUE, 20])
  minimumAL <- min(ALs)
  maximumAL <- max(ALs)
  
  ppy <- as.numeric(mapafit[1,21])
  
  if (fh == -1){
    fh <- ppy
  }
  
  # Override ifh.c to be >= fh if conf.lvl is not NULL
  # Output will crop values to ifh
  ifh.c <- ifh
  if (!is.null(conf.lvl)){
    if (ifh.c < fh){
      ifh.c <- fh
    } 
  }
  
  # In-sample MAPA
  i.start <- max(ppy,maximumAL)
  if (ifh.c>0 && i.start<observations){   # Do not produce in-sample forecasts if there is 
                                          # not enough sample
    infor <- array(NA,c(ifh.c,observations),dimnames=list(paste("t+",1:ifh.c,sep="")))
    for (i in i.start:(observations-1)){
      inobs <- as.matrix(y[1:i])
      infor[, i+1] <- mapacalc(inobs, mapafit, fh=ifh.c, comb, output="forecast", 
                               outplot=0, hybrid) 
      # Crop out-of-sample predictions
      if ((i+ifh.c)>observations){
        k <- (i+ifh.c) - observations
        infor[(ifh.c-k+1):ifh.c, i+1] <- rep(NA,k)
      }
    }    
  } else {
    infor <- NULL
    ifh.c <- 0        # Override in-sample output if no forecasts are calculated
    conf.lvl <- NULL  # Confidence intervals cannot be calculated
    
  }
  
  # Out-of-sample MAPA
  if (fh>0){
    outfor <- mapacalc(y, mapafit, fh, comb, output="forecast", outplot=0, hybrid) 
  } else {
    outfor <- NULL
  }
  
  # Calculate in-sample errors
  if (ifh.c == 1) {
    resid <- y - t(infor)
    MSE <- array(mean(resid^2, na.rm=TRUE),c(1,1),dimnames=list("t+1","MSE"))
    MAE <- array(mean(abs(resid), na.rm=TRUE),c(1,1),dimnames=list("t+1","MAE"))
  } else if (ifh.c > 1) {
    MSE <- array(NA,c(ifh.c,1),dimnames=list(paste("t+",1:ifh.c,sep=""),"MSE"))
    MAE <- array(NA,c(ifh.c,1),dimnames=list(paste("t+",1:ifh.c,sep=""),"MAE"))
    for (h in 1:min(ifh.c,(observations-ppy))) {
      resid <- y[h:observations] - infor[h, 1:(observations-h+1)]
      MSE[h] <- mean(resid^2, na.rm=TRUE)
      MAE[h] <- mean(abs(resid), na.rm=TRUE)
    }
  } else {
    MSE <- NULL
    MAE <- NULL
  }
  
  # Calculate prediction intervals
  if (!is.null(conf.lvl)){
    intv.idx <- !is.na(MSE)
    if (sum(!intv.idx)>0){
      intv <- c(sqrt(MSE[intv.idx]), sqrt(MSE[1])*sqrt((sum(intv.idx)+1):fh))
    } else {
      intv <- c(sqrt(MSE[1:fh]))
    }
    conf.lvl[conf.lvl>0.999999999999999] <- 0.999999999999999
    conf.lvl <- unique(conf.lvl)
    PIn <- length(conf.lvl)
    conf.lvl <- sort(conf.lvl,decreasing=TRUE)
    z <- abs(qnorm((1-conf.lvl)/2))
    PI <- array(NA,c(2*PIn,fh))
    for (i in 1:PIn){
      PI[i,] <- outfor + intv*z[i]
      PI[2*PIn-i+1,] <- outfor - intv*z[i]
    }
    rownames(PI) <- c(paste("Upper",format(conf.lvl,digit=2)),
                      paste("Lower",format(conf.lvl[PIn:1],digit=2)))
    colnames(PI) <- paste("t+",1:fh,sep="")
  } else {
    PI <- NULL
  }
  
  # Crop insample forecasts and erros to ifh
  if (ifh.c > 0 && ifh > 0){
    infor <- array(infor[1:ifh,],c(ifh,observations),dimnames=list(paste("t+",1:ifh,sep="")))
    MSE <- MSE[1:ifh,]
    MAE <- MAE[1:ifh,]
  } else {
    infor <- NULL
    MSE <- NULL
    MAE <- NULL
  }
  
  # Produce plot
  if (outplot==1){
    layout(matrix(1, 1, 1, byrow = TRUE))
    # Find min max
    if (is.null(outfor)){
      ymax <- max(y)
      ymin <- min(y)
      ymax <- ymax + 0.1*(ymax-ymin)
      ymin <- ymin - 0.1*(ymax-ymin)      
    } else {
      if (!is.null(conf.lvl)){
        ymax <- max(c(max(outfor),max(y),max(PI)))
        ymin <- min(c(min(outfor),min(y),min(PI)))
      } else {
        ymax <- max(c(max(outfor),max(y)))
        ymin <- min(c(min(outfor),min(y)))
      }
      ymax <- ymax + 0.1*(ymax-ymin)
      ymin <- ymin - 0.1*(ymax-ymin)
    }
    plot(1:observations,y,type="l",col="blue", xlab="", ylab="", main="Forecast", 
		xlim <- c(1, observations+fh), ylim=c(ymin,ymax))
    # In-sample
    if (ifh.c>0){
      if (ifh==1){
        lines(infor[1,],col="red")
      } else {
        # clrs = rainbow(observations-ppy)
        for (i in (ppy):(observations-1)){
          lines((i):(i+ifh-1),infor[,i],col="red")
        }
      }
    }
    # Prediction intervals
    if (!is.null(conf.lvl)){
      cmp <- rgb(1,0.75,0.75,1)
      for (i in 1:PIn){
        rr <- 1-(1-(i-1)*(0.5/PIn))*(1-c(col2rgb(cmp)/255))
        polygon(c(observations+(1:fh),observations+(fh:1)),c(PI[i,],PI[PIn*2+1-i,fh:1]),
                col=rgb(rr[1],rr[2],rr[3]), border=NA)
      }
    }
    # Out-of-sample
    if (ifh == 0 || ifh.c == 0){
      lines(observations:(fh+observations),c(y[observations],outfor),col="red")
    } else if (ifh == 1){
      lines(observations:(fh+observations),c(infor[1,observations],outfor),col="red")
    } else {
      lines((observations+1):(fh+observations),outfor,col="red")
    }
  }
  
  # Construct output
  return(list(infor=infor,outfor=outfor,PI=PI,MSE=MSE,MAE=MAE))
  
}

#-------------------------------------------------
mapaest <- function(y, ppy=NULL, minimumAL=1, maximumAL=ppy, paral=0, display=0, 
                    outplot=1, model="ZZZ") {
# Estimate MAPA for a time series  
#  
# Inputs:
#   y           = In sample observations of a time series (vector)
#   ppy         = Periods in a season of the time series at the sampled frequency.
#                 If y is a ts object then this is taken from its frequency,
#                 unless overriden. 
#   minimumAL   = Lowest aggregation level to use. Default = 1, maximumAL>1
#   maximumAL   = Highest aggregation level to use. Default = ppy
#   paral       = Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 
#                 2 = yes and initialise cluster. Default is 0.
#   display     = Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
#   outplot     = Provide output plot. 0 = no; 1 = yes. Default is 1.  
#   model       = Allow only that type of ETS at each aggregation level. This follows similar
#                 coding to the ets function. The first letter refers to the error 
#                 type ("A", "M" or "Z"); the second letter refers to the trend 
#                 type ("N","A","Ad","M","Md" or "Z"); and the third letter 
#                 refers to the season type ("N","A","M" or "Z"). The letters mean: 
#                 "N"=none, "A"=additive, "M"=multiplicative and "Z"=automatically selected. 
#                 A "d" for trend implies damped. By default model="ZZZ". If due to sample
#                 limitation ETS cannot be calculated at an aggregation level for the selected
#                 model, then no estimation is done for that specific level. For aggregation 
#                 levels that seasonality becomes 1 then a non-seasonal model is estimated.
#
# Output:
#   mapafit     = Estimated MAPA model structure

  # Get ppy and maximumAL
  if (is.null(ppy)){
    if (class(y)=="ts"){
      ppy <- frequency(y)
      if (is.null(maximumAL)){maximumAL <- ppy}
    } else {
      stop(paste("Input ppy is not given and y input not ts class.",
                 "Please provide the periods in a season of the time series",
                 "at the sampled frequency."))
    }
  }  
  
  # Make sure that maximumAL > 1 and larger than minimumAL
  if (maximumAL == 1){
    maximumAL = maximumAL + 1
  }
  
  if (minimumAL>=maximumAL){
    stop("maximumAL must be larger than minimumAL")
  }
  
  # Setup parallel processing if required
  if (paral == 2){
    crs <- detectCores()
    cl <- makeCluster(getOption("cl.cores", crs))
    writeLines(paste("Running with", crs, 'cores'))
  }
  
  observations <- length(y) # number of observations for the in-sample data
  
  # Aggregation and estimation
  if (paral != 0){  # Parallel run
    mapafit <- clusterApplyLB(cl, 1:(maximumAL-minimumAL+1), mapaest.loop, 
      y=y, minimumAL=minimumAL, maximumAL=maximumAL, observations=observations, ppy=ppy,
      display=display,model=model)  
  } else {          # Serial run
    mapafit <- vector("list", (maximumAL-minimumAL+1))
    for (i in 1:(maximumAL-minimumAL+1)){
      mapafit[[i]] <- mapaest.loop(i, y, minimumAL, maximumAL, observations, 
        ppy, display,model=model)
    }
  }
    
  if (paral == 2){
    # Stop parallel processing
    stopCluster(cl)
  }

  # Process output
  mapafit <- do.call(rbind, mapafit) # Re-arrange output for clusterApplyLB function
  rownames(mapafit) <- paste("AL",minimumAL:maximumAL,sep="")

  # mapafit <- mapafit[, c(20,19,11,12,13,14)]
  
  # Plot model selection summary
  ALplot <- 1:(maximumAL-minimumAL+1)
  ALplot <- ALplot[unlist(mapafit[,19])==TRUE]
  if (outplot == 1){
    layout(matrix(1, 1, 1, byrow = TRUE))
    comps <- array(0,c(max(ALplot),5))
    for (AL in 1:max(ALplot)){
      components <- mapafit[[AL, 14]]
      # Error term
      if (components[1]=="A"){
        comps[AL,1] <- 1
      } else {
        comps[AL,1] <- 2
      }
      # Trend term
      if (components[2]=="A"){
        comps[AL,2] <- 1
      } else {if (components[2]=="M"){
        comps[AL,2] <- 2
      } else
        comps[AL,2] <- 0
      }
      # Season term
      if (components[3]=="A"){
        comps[AL,3] <- 1
      } else {if (components[3]=="M"){
        comps[AL,3] <- 2
      } else
        comps[AL,3] <- 0
      }
      # Damped tem
      if (components[4]==TRUE){
        comps[AL,4] <- 1
      }
      comps[AL,5] <- mapafit[[AL,20]]
    }
    comps[, 2] <- comps[, 2] + 0.5*comps[, 4]
    image(min(comps[,5]):max(comps[,5]), 1:3, matrix(comps[,1:3],ncol=3), axes=FALSE, col=rev(heat.colors(5)), 
		  ylab="Components", xlab="Aggregation Level", main="ETS components")
    axis(2, at=1:3, labels=list("Error","Trend","Season"))
    axis(1, at=min(comps[,5]):max(comps[,5]))
    
    for (i in 1:4){
      for (AL in 1:max(ALplot)){
        if (i==1){
          lines(c(AL-0.5+minimumAL-1,AL-0.5+minimumAL-1),c(0,4),col="black")
        }
        if (i<4 & AL<=max(comps[,5])){
          if (i==2 & comps[AL,4]==TRUE){
            damp <- "d"
          } else {
            damp <- NULL
          }
          text(AL+minimumAL-1,i,paste(mapafit[[AL,14]][i],damp,sep=""))
        }
      }
      lines(c(min(comps[,5])-0.5,max(comps[,5])+0.5),c(i-0.5,i-0.5),col="black")
    }
    lines(c(as.numeric(mapafit[max(ALplot),20])+0.5,
            as.numeric(mapafit[max(ALplot),20])+0.5),c(0,4),col="black")
  }
  
  # Return output
  return(mapafit)
  
}

#-------------------------------------------------
mapacalc <- function(y, mapafit, fh=0, comb="mean", output="forecast", 
	outplot=0, hybrid=TRUE) 
{
# Calculation of MAPA forecasts
# 
# Inputs:
#   y           = In sample observations of a time series (vector)
#   mapafit     = Fitted MAPA model (from mapaest)
#   fh          = Forecast horizon. Default = ppy
#   comb        = Combination operator. One of "mean" or "median". Default is "mean"
#   output      = Type of output. One of "forecast" or "all". Default is "forecast"
#                 If output="all", both forecasts and components estimates per aggregation
#                 level are provided. For the components the rainbow colouring scheme is used. 
#                 Red is aggregation level 1, followed by yellow, green, cyan, blue and magenta 
#                 for the higher aggregation levels. 
#   outplot     = Provide output plot. 0 = no; 1 = time series and forecast only;
#                 2 = time series, forecasts and components. Default is 1. 
#   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE
#                 If minimumAL > 1 then the minimumAL ETS forecasts are used.
#
# Output:
#   forecasts   = Vector with forecasts
#   components  = array with MAPA components
  
  # Get settings from mapafit
  ALs <- as.numeric(mapafit[mapafit[,19]==TRUE, 20])
  minimumAL <- min(ALs)
  maximumAL <- max(ALs)
  ppy <- as.numeric(mapafit[1,21])
  
  # Set default foreast horizon
  if (fh == 0){
    fh <- ppy
  }
    
  observations <- length(y) # number of observations for the in-sample data
  
  FCs <- array(0, c(maximumAL-minimumAL+1, 4, fh),dimnames=list(paste("AL",minimumAL:maximumAL,sep=""),
	c("ETS","Level","Trend","Season"),paste("t+",1:fh,sep=""))) # the forecasted components 
																# are saved here
  
  # MAPA forecast
  ALvec <- minimumAL:maximumAL
  
  for (ALi in 1:(maximumAL-minimumAL+1)){
    
    AL <- ALvec[ALi]
    
    q <- observations %/% AL # observation in the aggregated level
    r <- observations %% AL  # observation to discard from the beginning of the series
    ppyA <- ppy %/% AL       # periods per year for the aggregated level
    if (ppy %% AL != 0){
      ppyA <- 1
    }
    
    # Aggregation
    yA <- array(0, dim=c(q)) # in-sample aggregated values will be saved here
    for (j in 1:q){                 # calculate the aggregate values
      yA[j] <- mean(y[(r+1+(j-1)*AL):(r+j*AL)])
    }
    
    ats <- ts(yA, frequency = ppyA) # create the time series

# This part of code used to call forecast:::pegelsresid.C
# This is now obsolete with the new forecast package that has use.initial.values=TRUE
#
#     # Extarct ets fit from mapafit
#     components <- mapafit[[AL,14]]
#     errortype <- components[1]
#     trendtype <- components[2]
#     seasontype <- components[3]
#     damped <- as.logical(components[4])
#     param <- mapafit[[AL,11]]
#     pnames <- names(param)
#     alpha <- param[pnames=="alpha"]
#     if (sum(pnames=="beta")>0) {
#       beta <- param[pnames=="beta"]
#     } else {
#       beta <- NULL
#     }
#     if (sum(pnames=="gamma")>0) {
#       gamma <- param[pnames=="gamma"]
#     } else {
#       gamma <- NULL
#     }
#     if (sum(pnames=="phi")>0) {
#       phi <- param[pnames=="phi"]
#     } else {
#       phi <- NULL
#     }
#     
#     # Prepare initial states
#     np <- length(param)
#     nstate <- np - length(c(alpha,beta,gamma,phi))
#     initstate <- param[(np-nstate+1):np]
#     # Add extra state as per ets.R
#     if(seasontype!="N")
#       initstate <- c(initstate, ppyA*(seasontype=="M")-sum(initstate[(2+(trendtype!="N")):nstate]))
#     
#     # Calculate ets states using mapafit results
#     ats.fit <- forecast:::pegelsresid.C(ats, ppyA, initstate, errortype, trendtype,
# 		  seasontype, damped, alpha, beta, gamma, phi)
#     ats.fit$components <- components
    
    # ETS based calculation
    param <- mapafit[[ALi,11]]
    pnames <- names(param)
    if (sum(pnames=="phi")>0) {
      phi <- param[pnames=="phi"]
    } else {
      phi <- NULL
    }
    
    AL.fit <- structure(mapafit[ALi,1:18],class="ets")
    ats.fit <- ets(ats, AL.fit, use.initial.values=TRUE)
    
    # Transalte ets states for MAPA
    FCs_temp <- statetranslate(ats.fit,AL,fh,q,ppyA,phi,1)
    
    # Return MAPA components
    FCs[ALi, , ] <- FCs_temp
    
  }
  
  # MAPA combination
  combres <- mapacomb(minimumAL,maximumAL,ppy,FCs,comb,observations)
  forecasts <- combres$forecasts
  perm_levels <- combres$perm_levels
  perm_seas <- combres$perm_seas

  # Calculate hybrid model
  if (hybrid==TRUE){
    forecasts <- (FCs[1,1,] + forecasts)/2
  }  
  
  # Plot output
  mapaplot(outplot,FCs,maximumAL,perm_levels,perm_seas,observations,y,forecasts,fh,comb)

  # Construct output
  if (output=="forecast"){
    return(forecasts)
  } else {
    return(list(forecast=forecasts,components=FCs))
  }
  
}

#-------------------------------------------------
statetranslate <- function(fit,AL,fh,q,ppyA,phi,fittype){
# This function prepares ets states for MAPA combination
# It extrapolates from last states the forecasts and translates to additive
  
  FCs_temp <- array(0, c(4, fh))
  
  fhA <- (fh %/% AL) + 1   # forecast horizon for the aggregated level
  
  # Estimates for the Level Component
  FCs_temp[2, ] <- as.numeric(rep(rep(fit$states[q+1, 1], fhA), each=AL)[1:fh])
  
  # Estimates for the Trend Component
  if (fit$components[2]=="N"){ # no trend
    FCs_temp[3, ] <- 0
    b = 0 # indicates that there is no trend
  } else if (fit$components[2]=="A"){ # additive trend
    if (fit$components[4]=="FALSE"){ 
      FCs_temp[3, ] <- as.numeric(rep(fit$states[q+1, 2] * (1:fhA), each=AL))[1:fh]
    } else { # additive damped trend
      # We divide with phi because of an internal calculation for 
	  # the damped trend in the ETS package
      FCs_temp[3, ] <- as.numeric(rep(cumsum((fit$states[q+1, 2]/phi)* phi^(1:fhA)), each=AL))[1:fh]
    }
    b <- 1 # indicates that there is trend
  } else {
    if (fit$components[4]=="FALSE"){ # multiplicative trend
      FCs_temp[3, ] <- as.numeric(rep((fit$states[q+1,2]^(1:fhA)-1), each=AL)[1:fh] * FCs_temp[2,])
    } else { # multiplicative damped trend
      # We divide with phi because of an internal calculation for the damped trend in the ETS package
      FCs_temp[3, ] <- as.numeric(rep((((fit$states[q+1,2] ^ (1/phi)) ^ cumsum(phi^(1:fhA)))-1), 
		each=AL)[1:fh] * FCs_temp[2, ])
    }
    b <- 1 # indicates that there is trend
  }
  
  # Estimates for the Seasonal Component 
  if (fit$components[3]=="N"){ # no seasonality
    FCs_temp[4, ] <- 0
  } else if (fit$components[3]=="A"){ # additive seasonality
    FCs_temp[4, ] <- as.numeric(rep(rep(rev(fit$states[q+1,(2+b):(ppyA+1+b)]), fhA), each=AL))[1:fh]
  } else { # multiplicative seasonality
    FCs_temp[4, ] <- as.numeric((rep(rep(rev(fit$states[q+1,(2+b):(ppyA+1+b)]), fhA),
		                 each=AL)[1:fh] - 1)) * (FCs_temp[2, ] + FCs_temp[3, ])
  }
  
  # fittype identifies if information is comming from ets or mapafit
  if (fittype==1){
    # Recreate ETS forecasts
    if (fh != 1) {
      FCs_temp[1, ] <- colSums(FCs_temp[2:4,])     
    } else {
      FCs_temp[1, ] <- sum(FCs_temp[2:4,])   
    }
  }

  # Return output
  return(FCs_temp)
  
}

#-------------------------------------------------
mapacomb <- function(minimumAL,maximumAL,ppy,FCs,comb,observations){
# This function combines the translated ets states
  
  # perm_levels is not needed for forecasting. This is already checked in the estimation.
  # perm_levels <- array(0, maximumAL) # permitted levels due to ETS implementation (observations>=4)
  perm_seas <- array(0, maximumAL)   # permitted seasonalities
  for (AL in minimumAL:maximumAL){
    # if (observations %/% AL >=4){
    #   perm_levels[AL] <- 1
    # }
    if ((ppy %% AL == 0) & (AL<ppy)){
      perm_seas[AL] <- 1
    }
  }
  # perm_levels <- perm_levels[minimumAL:maximumAL]
  perm_levels <- rep(1,(maximumAL-minimumAL+1))
  perm_seas <- perm_seas[minimumAL:maximumAL]
  
  if (dim(FCs)[3] != 1){ # Forecast multiple steps ahead
    level <- FCs[perm_levels==1, 2, ]
    trend <- FCs[perm_levels==1, 3, ]
    season <- FCs[(perm_levels==1 & perm_seas==1), 4, ]
    # Check that all are arrays
    if (!is.array(level)){
      level <- array(level,c(1,length(level)))
    }
    if (!is.array(trend)){
      trend <- array(trend,c(1,length(trend)))
    }
    if (!is.array(season)){
      season <- array(season,c(1,length(season)))
    }
    if (comb=="mean"){ # alternative averaging operators
      forecasts <- colSums(rbind(colMeans(level),colMeans(trend),
                                 colMeans(season)),na.rm=TRUE) # MAPA(mean) forecasts
    } else {
      # forecasts <- colSums(rbind(colMedians(level),colMedians(trend),
      #                            colMedians(season)),na.rm=TRUE) # MAPA(median) forecasts
      forecasts <- colSums(rbind(apply(level,2,"median"),apply(trend,2,"median"),
                                 apply(season,2,"median")),na.rm=TRUE) # MAPA(median) forecasts
    }
  } else {
    if (comb=="mean"){ # alternative averaging operators
      forecasts <- colSums(rbind(mean(FCs[perm_levels==1, 2, ]),mean(FCs[perm_levels==1, 3, ]),
		    mean(FCs[(perm_levels==1 & perm_seas==1), 4, ])), na.rm=TRUE) # MAPA(mean) forecasts
    } else {
      forecasts <- colSums(rbind(median(FCs[perm_levels==1, 2, ]),median(FCs[perm_levels==1, 3, ]),
		    median(FCs[(perm_levels==1 & perm_seas==1), 4, ])), na.rm=TRUE) # MAPA(mean) forecasts
    }
  }
  
  # Return output
  return(list(forecasts=forecasts,perm_levels=perm_levels,perm_seas=perm_seas))
}

#-------------------------------------------------
mapaplot <- function(outplot,FCs,maximumAL,perm_levels,perm_seas,observations,
	                   y,forecasts,fh,comb)
{
# Produce MAPA forecast & components plot 
# outplot == 0, no plot; == 1 series plot; == 2 component plot  

  if (outplot > 0){
    FClevel <- FCs[perm_levels==1, 2, ]
    FCtrend <- FCs[perm_levels==1, 3, ]
    if (sum(perm_levels==1 & perm_seas==1)!=0){
      FCseason <- FCs[(perm_levels==1 & perm_seas==1), 4, ]
    } else {
      FCseason <- NULL
    }
    # Check that all are arrays
    if (!is.array(FClevel)){
      FClevel <- array(FClevel,c(length(FClevel)/fh,fh))
    }
    if (!is.array(FCtrend)){
      FCtrend <- array(FCtrend,c(length(FCtrend)/fh,fh))
    }
    if (!is.array(FCseason) && !is.null(FCseason)){
      FCseason <- array(FCseason,c(length(FCseason)/fh,fh))
    }
    clrs <- rainbow(length(perm_levels))
    if (outplot == 2){
      if (!is.null(FCseason)){
        layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))
      } else {
        layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
      }
    } else {
      layout(matrix(1, 1, 1, byrow = TRUE))
    }
    # Find min max
    ymax <- max(c(max(forecasts),max(y)))
    ymin <- min(c(min(forecasts),min(y)))
    yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
    # Plot prediction
    plot(1:observations,y, type="l", col="blue", xlab="", ylab="", 
		  main="Forecast", xlim=c(1, observations+fh), ylim=yminmax)
    lines((observations):(observations+fh),c(y[observations],forecasts), col="red")
    if (outplot == 2){
      # Plot Level
      ymin <- min(FClevel)
      ymax <- max(FClevel)
      yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
      plot(FClevel[1, ], type="l", col=clrs[1], xlab="", ylab="", main="Level", ylim=yminmax)
      if (sum(perm_levels)>1){
        for (i in 2:sum(perm_levels)){
          lines(FClevel[i, ], type="l", col=clrs[i])
        }
      }
      if (comb=="mean"){
        lines(colMeans(FClevel), type="l", col="black", lwd=2)
      } else {
        # lines(colMedians(FClevel), type="l", col="black", lwd=2)
        lines(apply(FClevel,2,"median"), type="l", col="black", lwd=2)
      }
      # Plot trend
      ymin <- min(FCtrend)
      ymax <- max(FCtrend)
      yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
      plot(FCtrend[1, ], type="l", col=clrs[1], xlab="", ylab="", main="Trend", 
		    ylim=yminmax)
      if (sum(perm_levels)>1){
        for (i in 2:sum(perm_levels)){
          lines(FCtrend[i, ], type="l", col=clrs[i])
        }
      }
      if (comb=="mean"){
        lines(colMeans(FCtrend), type="l", col="black", lwd=2)
      } else {
        # lines(colMedians(FCtrend), type="l", col="black", lwd=2)
        lines(apply(FCtrend,2,"median"), type="l", col="black", lwd=2)
      }
      # Plot season
      if (!is.null(FCseason)){
        ymin <- min(FCseason)
        ymax <- max(FCseason)
        yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
        plot(FCseason[1, ], type="l", col=clrs[1], xlab="", ylab="", main="Season", ylim=yminmax)
        if (dim(FCseason)[1]>1){
          for (i in 2:sum(perm_levels==1 & perm_seas==1)){
            lines(FCseason[i, ], type="l", col=clrs[i])
          }
        }
        if (comb=="mean"){
          lines(colMeans(FCseason), type="l", col="black", lwd=2)
        } else {
          # lines(colMedians(FCseason), type="l", col="black", lwd=2)
          lines(apply(FCseason,2,"median"), type="l", col="black", lwd=2)
        }
      }
    }
  }
}

#-------------------------------------------------
mapasimple.loop <- function(ALi, y, minimumAL, maximumAL, observations, 
                            ppy, display, fh, model){
  # Internal function for mapasimple estimation and forecast iterations
  
  # Create ETS model strings
  mn <- nchar(model)
  model <- substring(model, seq(1,mn,1), seq(1,mn,1))
  ets.model <- paste(model[1],model[2],model[mn],sep="")
  if (model[2]=="Z"){
    ets.damped <- NULL
  } else if (mn==4){
    ets.damped <- TRUE
  } else {
    ets.damped <- FALSE
  }
  
  ALvec <- minimumAL:maximumAL
  AL <- ALvec[ALi]
  
  # Console output
  if (display==1){
    txtc <- paste("Aggregation level: ",AL,"/",maximumAL,
                  " (",round(100*ALi/(maximumAL-minimumAL+1),2),"%)",sep="")
    cat(txtc)
  }
  
  FCs_temp <- array(NA, c(4, fh))
  
  q <- observations %/% AL # observation in the aggregated level
  r <- observations %% AL  # observation to discard from the beginning of the series
  fhA <- (fh %/% AL) + 1   # forecast horizon for the aggregated level
  ppyA <- ppy %/% AL       # periods per year for the aggregated level
  if (ppy %% AL != 0){
    ppyA <- 1
  }
  
  # Check if selected ETS model is possible for current AL
  npars <- 2
  if (model[2] == "A" | model[2] == "M"){
    npars <- npars + 2}
  if (model[mn] == "A" | model[mn] == "M"){
    npars <- npars + ppyA}
  if (!is.null(ets.damped)){ 
    npars <- npars + as.numeric(ets.damped)}
  if (q <= npars + 1){
    q <- 1} # This will not estimate current AL
  
  if (q >= 4){
    # Aggregation
    yA <- array(0, dim=c(q)) # in-sample aggregated values will be saved here
    for (j in 1:q){                 # calculate the aggregate values
      yA[j] <- mean(y[(r+1+(j-1)*AL):(r+j*AL)])
    }
    
    ats <- ts(yA, frequency = ppyA) # create the time series
    
    # Check if seasonality exists and select appropriate model
    if ((model[mn] == "A" | model[mn] == "M") && ppyA == 1){
      mapa.model <- paste(substring(ets.model,1,2),"N",sep="")
    } else {
      mapa.model <- ets.model
    }
    
    # Fit ETS
    ats.fit <- ets(ats,model=mapa.model,damped=ets.damped)
    ats.fcast <- forecast.ets(ats.fit,h=fhA)
    
    # Translate ets states for MAPA
    phi <- ats.fit$par[names(ats.fit$par)=="phi"]
    FCs_temp <- statetranslate(ats.fit,AL,fh,q,ppyA,phi,0)
    
    # ets forecast on original frequency
    FCs_temp[1, ] <- rep(ats.fcast$mean[1:fhA], each=AL)[1:fh]
  }
  
  # Update console display
  if (display==1){
    nc <- nchar(txtc)
    cat(rep("\r",nc))
    cat(rep(" ",nc))
    cat(rep("\r",nc))
  }
  
  # Return FCs_temp for foreach
  return(FCs_temp)
  
}

#-------------------------------------------------
mapaest.loop <- function(ALi, y, minimumAL, maximumAL, observations, 
                         ppy, display, model){ 
  # Internal function for running a single loop in mapaest
  
  # Create ETS model strings
  mn <- nchar(model)
  model <- substring(model, seq(1,mn,1), seq(1,mn,1))
  ets.model <- paste(model[1],model[2],model[mn],sep="")
  if (model[2]=="Z"){
    ets.damped <- NULL
  } else if (mn==4){
    ets.damped <- TRUE
  } else {
    ets.damped <- FALSE
  }
  
  ALvec <- minimumAL:maximumAL
  AL <- ALvec[ALi]
  
  # Console output
  if (display==1){
    txtc <- paste("Aggregation level: ",AL,"/",maximumAL,
                  " (",round(100*ALi/(maximumAL-minimumAL+1),2),"%)",sep="")
    cat(txtc)
  }
  
  q <- observations %/% AL # observation in the aggregated level
  r <- observations %% AL  # observation to discard from the beginning of the series
  ppyA <- ppy %/% AL       # periods per year for the aggregated level
  if (ppy %% AL != 0){
    ppyA <- 1
  }
  
  # Check if selected ETS model is possible for current AL
  npars <- 2
  if (model[2] == "A" | model[2] == "M"){
    npars <- npars + 2}
  if (model[mn] == "A" | model[mn] == "M"){
    npars <- npars + ppyA}
  if (!is.null(ets.damped)){ 
    npars <- npars + as.numeric(ets.damped)}
  if (q <= npars + 1){
    q <- 1} # This will not estimate current AL
  
  if (q >= 4){
    # Aggregation
    yA <- array(0, dim=c(q)) # in-sample aggregated values will be saved here
    for (j in 1:q){                 # calculate the aggregate values
      yA[j] <- mean(y[(r+1+(j-1)*AL):(r+j*AL)])
    }
    
    ats <- ts(yA, frequency = ppyA) # create the time series
    
    # Check if seasonality exists and select appropriate model
    if ((model[mn] == "A" | model[mn] == "M") && ppyA == 1){
      mapa.model <- paste(substring(ets.model,1,2),"N",sep="")
    } else {
      mapa.model <- ets.model
    }
    
    # Fit ETS
    fit <- ets(ats,model=mapa.model,damped=ets.damped)
    # If time series is constant then ets does not return the same
    # structure, correct for that
    if (length(names(fit)) != 18){
      fit.temp = fit
      fit <- list("loglik"=NULL,"aic"=NULL,"bic"=NULL,"aicc"=NULL,"mse"=NULL,
                  "amse"=NULL,"fit"=NULL,"residuals"=NULL,"fitted"=NULL,
                  "states"=NULL,"par"=NULL,"m"=NULL,"method"=NULL,
                  "components"=NULL,"call"=NULL,"initstate"=NULL,
                  "sigma2"=NULL,"x"=NULL)
      fit.temp.names <- names(fit.temp)
      fit.names <- names(fit)
      for (fi in 1:length(fit.names)){
        if (sum(fit.temp.names == fit.names[fi])>0){
          eval(parse(text=paste0("fit$",fit.names[fi]," <- fit.temp$",fit.names[fi])))
        }
      }
    }
    fit$use <- TRUE
  } else {
    fit <- list("loglik"=NULL,"aic"=NULL,"bic"=NULL,"aicc"=NULL,"mse"=NULL,
                "amse"=NULL,"fit"=NULL,"residuals"=NULL,"fitted"=NULL,
                "states"=NULL,"par"=NULL,"m"=NULL,"method"=NULL,
                "components"=NULL,"call"=NULL,"initstate"=NULL,
                "sigma2"=NULL,"x"=NULL,"use"=FALSE)
  }
  
  fit$AL <- AL
  fit$original.ppy <- ppy
  
  # Update console display
  if (display==1){
    nc = nchar(txtc)
    cat(rep("\r",nc))
    cat(rep(" ",nc))
    cat(rep("\r",nc))
  }
  
  # Return loop result
  return(rbind(fit))
  
}

#-------------------------------------------------
plotmapa <- function(mapafit){
# Produce estimated MAPA fit plot
# 
# Inputs:
#   mapafit     = Fitted MAPA model (from mapaest)


  # Get settings from mapafit
  ALs <- as.numeric(mapafit[mapafit[,19]==TRUE, 20])
  minimumAL <- min(ALs)
  maximumAL <- max(ALs)
  ppy <- as.numeric(mapafit[1,21])
  
  # Plot model selection summary
  ALplot <- 1:(maximumAL-minimumAL+1)
  ALplot <- ALplot[unlist(mapafit[,19])==TRUE]

  layout(matrix(1, 1, 1, byrow = TRUE))
  comps <- array(0,c(max(ALplot),5))
  for (AL in 1:max(ALplot)){
    components <- mapafit[[AL, 14]]
    # Error term
    if (components[1]=="A"){
      comps[AL,1] <- 1
    } else {
      comps[AL,1] <- 2
    }
    # Trend term
    if (components[2]=="A"){
      comps[AL,2] <- 1
    } else {if (components[2]=="M"){
      comps[AL,2] <- 2
    } else
      comps[AL,2] <- 0
    }
    # Season term
    if (components[3]=="A"){
      comps[AL,3] <- 1
    } else {if (components[3]=="M"){
      comps[AL,3] <- 2
    } else
      comps[AL,3] <- 0
    }
    # Damped tem
    if (components[4]==TRUE){
      comps[AL,4] <- 1
    }
    comps[AL,5] <- mapafit[[AL,20]]
  }
  comps[, 2] <- comps[, 2] + 0.5*comps[, 4]
  image(min(comps[,5]):max(comps[,5]), 1:3, matrix(comps[,1:3],ncol=3), axes=FALSE, col=rev(heat.colors(5)), 
        ylab="Components", xlab="Aggregation Level", main="ETS components")
  axis(2, at=1:3, labels=list("Error","Trend","Season"))
  axis(1, at=min(comps[,5]):max(comps[,5]))
  
  for (i in 1:4){
    for (AL in 1:max(ALplot)){
      if (i==1){
        lines(c(AL-0.5+minimumAL-1,AL-0.5+minimumAL-1),c(0,4),col="black")
      }
      if (i<4 & AL<=max(comps[,5])){
        if (i==2 & comps[AL,4]==TRUE){
          damp <- "d"
        } else {
          damp <- NULL
        }
        text(AL+minimumAL-1,i,paste(mapafit[[AL,14]][i],damp,sep=""))
      }
    }
    lines(c(min(comps[,5])-0.5,max(comps[,5])+0.5),c(i-0.5,i-0.5),col="black")
  }
  lines(c(as.numeric(mapafit[max(ALplot),20])+0.5,
          as.numeric(mapafit[max(ALplot),20])+0.5),c(0,4),col="black")
  
}