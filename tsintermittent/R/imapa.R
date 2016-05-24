# Based on paper:
# Petropoulos F. and Kourentzes N. (2014) 
# Forecast Combinations for Intermittent Demand.
# Journal of Operational Research Society.

library(MAPA, quietly=TRUE) 	   # Needed for tsaggr
library(parallel, quietly=TRUE)  # Needed for parallel

#-------------------------------------------------

imapa <- function(data,h=10,w=NULL,minimumAL=1,maximumAL=NULL,comb=c("mean","median"),
                  init.opt=c(TRUE,FALSE),paral=c(0,1,2),outplot=c(0,1,2,3,4),
                  model.fit=NULL,na.rm=c(FALSE,TRUE)){
# MAPA for intermittent demand data
#
# Inputs
#   data        Intermittent demand time series.
#   h           Forecast horizon.
#   w           Smoothing parameters. If w == NULL then parameters are optimised.
#               If w is w single parameter then the same is used for smoothing both the 
#               demand and the intervals. If two parameters are provided then the second 
#               is used to smooth the intervals. SES is always optimised.
#   minimumAL   Lowest aggregation level to use. Default = 1
#   maximumAL   Highest aggregation level to use. Default = maximum interval
#   comb        Combination operator. One of "mean" or "median". Default is "mean"
#   init.opt    If init.opt==TRUE then Croston and SBA initial values are optimised. 
#   paral       Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 
#               2 = yes and initialise cluster. Default is 0.
#   outplot     Optional plot:
#                 0 - No plot
#                 1 - Time series and combined forecast
#                 2 - Time series and all aggregation level forecasts 
#                 3 - Summary model selection plot
#                 4 - Detailed model selection plot
#   model.fit   Optional input with model types and parameters. This is the model.fit 
#               output from this function. If used it overrides other model settings.
#   na.rm       A logical value indicating whether NA values should be remove using the method.    
#
# Outputs
#   frc.in      In-sample demand rate. 
#   frc.out     Out-of-sample demand rate.
#   summary     An array containing information for each aggregation level:
#                 AL - Aggregation level
#                 n - Number of observations of aggregated series
#                 p - Average inter-demand interval
#                 cv2 - Coefficient of variation squared of non-zero demand
#                 model - Selected model, where 1 is Croston, 2 is SBA and 3 is SES
#                 use - If == 0 then this aggregation level is ignored because it 
#                       contains less than 4 observations.
#   model.fit   Parameters and initialisation values of fitted model in each aggregation
#               level, with aggregation level and model information.
#
# Example:
#   imapa(ts.data1,outplot=1)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>

  comb <- comb[1]
  paral <- paral[1]
  outplot <- outplot[1]
  init.opt <- init.opt[1]
  na.rm <- na.rm[1]

  # Prepare data
  if (class(data)=="data.frame"){
    if (ncol(data)>1){
      warning("Data frame with more than one columns. Using only first one.")
    }
    data <- data[[1]]
  }
  if (na.rm == TRUE){
    data <- data[!is.na(data)]
  }
  n <- length(data)
  
  # Check number of non-zero values - need to have at least two
  if (sum(data!=0)<2){
    stop("Need at least two non-zero values to model time series.")
  }
  
  # Setup parallel processing if required
  if (paral == 2){
    crs <- detectCores()
    cl <- makeCluster(getOption("cl.cores", crs))
    writeLines(paste("Running with", crs, 'cores'))
    invisible(clusterCall(cl, function(pkgs) {
      require(tsintermittent)
    }))
  }
  
  # Find interval alpha
  if (!is.null(w)){
    w.in <- w[length(w)]
  } else {
    w.in <- NULL
  }
  
  # If model.fit is provided override minimum and maximum AL
  if (!is.null(model.fit)){
    minimumAL <- min(model.fit[1,])
    maximumAL <- max(model.fit[1,])
  } 
  
  # If no maximumAL find maximum interval
  if (is.null(maximumAL)){
    maximumAL <- floor(n/2)
  }
  
  # Check minimumAL and maximumAL
  if (minimumAL>maximumAL){
    stop("minimumAL must be smaller or equal to maximumAL.")
  }
  
  # Aggregate time series
  yaggr <- tsaggr(data,minimumAL:maximumAL)
  
  # Model each aggregation level
  k <- length(yaggr$out)
  ychar <- array(NA,c(6,k),dimnames=list(c("AL","n","p","cv2","model","use"), names(yaggr$out)))
  fit <- array(NA,c(6,k),dimnames=list(c("AL","model","a1","a2","initial.z","initial.x"), names(yaggr$out)))
  f.in <- array(NA,c(k,n),dimnames=list(names(yaggr$out)))
  f.out <- array(NA,c(k,h),dimnames=list(names(yaggr$out)))
  
  # Forecast calculation
  if (paral != 0){  # Parallel run
    Fc_par <- clusterApplyLB(cl, 1:k, imapa.loop, yaggr=yaggr, minimumAL=minimumAL, w=w, 
                             w.in=w.in, init.opt, n=n, h=h, model.fit=model.fit)
  } else {          # Serial run
    Fc_par <- vector("list", k)
    for (i in 1:k){
      Fc_par[[i]] <- imapa.loop(i,yaggr,minimumAL,w,w.in,init.opt,n,h,model.fit)
    }
  }
  
  # Distribute parallel output
  Fc_par <- do.call(rbind, Fc_par)
  ychar[] <- t(Fc_par[,1:6])
  fit[] <- t(Fc_par[,7:12])
  f.out[] <- Fc_par[,13:(12+h)]
  f.in[] <- Fc_par[,(13+h):length(Fc_par[1,])]
  
  if (paral == 2){
    # Stop parallel processing
    stopCluster(cl)
  }
  
  # Combine forecasts
  if (sum(ychar[6,])>1){ 
    # If multiple aggregation levels are used
    if (comb=="median"){ # median
      frc.in <- apply(f.in[ychar[6,]==1,],2,"median")
      frc.out <- apply(f.out[ychar[6,]==1,],2,"median")
    } else { # mean
      frc.in <- colMeans(f.in[ychar[6,]==1,])
      frc.out <- colMeans(f.out[ychar[6,]==1,])
    }
  } else {
    # Single aggregation level
    frc.in <- f.in # f.in[ychar[6,]==1,]
    frc.out <- f.out # f.out[ychar[6,]==1,]
  }
  
  # Produce plots
  if (outplot==1){
    # Summary forecast
    plot(1:n,as.vector(data),type="l",xlim=c(1,(n+h)),xlab="Period",ylab="",
         xaxs="i",yaxs="i",ylim=c(0,max(data)*1.1))
    lines(which(data>0),data[data>0],type="p",pch=20)
    lines(1:n,frc.in,col="red")
    lines((n+1):(n+h),frc.out,col="red",lwd=2)
  } 
  if (outplot==2){
    # Detail forecast
    cmp = rainbow(sum(ychar[6,]),start=2/6,end=4/6)
    plot(1:n,as.vector(data),type="l",xlim=c(1,(n+h)),xlab="Period",ylab="",
         xaxs="i",yaxs="i",ylim=c(0,max(data)*1.1),lwd=2)
    lines(which(data>0),data[data>0],type="p",pch=20,lwd=2)
    for (i in 1:sum(ychar[6,])){
      lines(1:n,f.in[i,],col=cmp[i])
      lines((n+1):(n+h),f.out[i,],col=cmp[i])
    }
    lines(1:n,frc.in,col="red",lwd=2)
    lines((n+1):(n+h),frc.out,col="red",lwd=2)
    if (k>1){
      legend("topleft",c(paste0("AL",ychar[1,1]),paste0("AL",ychar[1,i]),"Combined"),
             col=c(cmp[1],cmp[i],"red"),lty=1,lwd=c(1,1,2),cex=0.6)
    } else {
      legend("topleft",c(paste0("AL",ychar[1,1]),"Combined"),
             col=c(cmp[1],"red"),lty=1,lwd=c(1,2),cex=0.6)
    }
  }
  if (outplot==3){
    # Model selection summary
    xx <- c(0.94, 1.5)
    yy <- c(0, 1)
    plot(0, 0, xlim = xx, ylim = yy, xaxs = "i", yaxs = "i", 
         xlab = "p", ylab = parse(text = "CV^2"), xaxt = "n", yaxt = "n")
    p.x <- c(seq(1, 4/3, 0.02), 4/3)
    w.p <- 1
    v.data <- (4*p.x*(2-p.x)-w.p*(4-w.p)-p.x*(p.x-1)*(4-w.p)*(2-w.p))/(p.x*(4-w.p)*(2*p.x-w.p))
    polygon(c(p.x, 1), c(v.data, 0), col = "lightgrey", border = NA)
    w.p <- 0
    v.y2 <- (4*p.x*(2-p.x)-w.p*(4-w.p)-p.x*(p.x-1)*(4-w.p)*(2-w.p))/(p.x*(4-w.p)*(2*p.x-w.p))
    v.y2[v.y2 < 0] <- 0
    polygon(c(p.x, p.x[18:1]), c(v.data, v.y2[18:1]), col = gray(0.2,0.3), border = NA)
    text(xx[2], 1, paste("SBA: ", sum(ychar[5,]==2,na.rm=TRUE), sep = ""),adj = c(1.2, 1.2))
    text(1, 0, paste("Croston: ", sum(ychar[5,]==1,na.rm=TRUE), sep = ""),adj = c(-0.2, -1.2))
    polygon(c(0, 1, 1, 0), c(0, 0, yy[2], yy[2]), col = gray(0.6), border = NA)
    text(0.97, 0, paste("SES: ", sum(ychar[5,]==3,na.rm=TRUE), sep = ""), adj = c(-0.2, 0.2), srt = 90)
  }
  if (outplot==4){
    # Model selection detail
    kmax <- max(which(ychar[6,]==1))
    cmp = rainbow(kmax,start=2/6,end=4/6)
    p <- ychar[3,ychar[6,]==1]
    v <- ychar[4,ychar[6,]==1]
    xmin <- 1 - (max(c(max(p) + diff(range(p)) * 0.1), 1.5) - 1) * 0.1
    xx <- c(xmin, max(c(max(p) + diff(range(p)) * 0.1), 1.5))
    yy <- c(0, max(c(max(v) + diff(range(v)) * 0.1), 1.5))
    plot(0, 0, type="p", pch=20, xaxs="i", yaxs="i", xlim=xx, ylim=yy, xlab="p", 
         ylab=parse(text="CV^2"))
    p.x <- c(seq(1, 4/3, 0.02), 4/3)
    w.p <- 1
    v.data <- (4*p.x*(2-p.x)-w.p*(4-w.p)-p.x*(p.x-1)*(4-w.p)*(2-w.p))/(p.x*(4-w.p)*(2*p.x-w.p))
    polygon(c(p.x, 1), c(v.data, 0), col = "lightgrey", border = NA)
    w.p <- 0
    v.y2 <- (4*p.x*(2-p.x)-w.p*(4-w.p)-p.x*(p.x-1)*(4-w.p)*(2-w.p))/(p.x*(4-w.p)*(2*p.x-w.p))
    v.y2[v.y2 < 0] <- 0
    polygon(c(p.x, p.x[18:1]), c(v.data, v.y2[18:1]), col = gray(0.2,0.3), border = NA)
    polygon(c(0, 1, 1, 0), c(0, 0, yy[2], yy[2]), col = gray(0.4), 
            border = NA)
    
    lines(p[p>1], v[p>1], type="p", pch=20, col=cmp[p>1])
    lines(rep(xmin+(1-xmin)/2,sum(p==1)), v[p==1], type="p", pch=20, col=cmp[p==1])
    if (k>1){
      legend("topright",c(paste0("AL",ychar[1,1]),paste0("AL",ychar[1,kmax])),
             col=c(cmp[1],cmp[kmax]),pch=20,cex=0.6)
    } else {
      legend("topright",paste0("AL",ychar[1,1]),col=cmp[1],pch=20,cex=0.6)
    }
  }
 
  return(list(frc.in=frc.in,frc.out=frc.out,summary=ychar,model.fit=fit))
  
}

#-------------------------------------------------

imapa.loop <- function(i,yaggr,minimumAL,w,w.in,init.opt,n,h,model.fit){
# Forecast for each individual aggregation level  - internal function
  
  ytemp <- as.vector(yaggr$out[[i]])
  ychartemp <- array(NA,c(1,6))
  
  # Get time series characteristics
  ychartemp[1] <- minimumAL + i - 1
  ychartemp[2] <- length(ytemp)
  ychartemp[6] <- ychartemp[2]>=4
  
  # Check if at this aggregation level there are at least two non-zero observations
  if (ychartemp[6]==1 && sum(ytemp!=0)<2){
    ychartemp[6] <- 0
  }
  
  if (ychartemp[6]==1){
    if (!is.null(model.fit)){
      # Select prefit model
      fit <- model.fit[,i]
      if (!is.na(fit[2])){
        # chartemp <- idclass(ytemp,type="PK",a.in=0.1,outplot="none")
        chartemp <- idclass(ytemp,type="PKa",a.in=0.1,outplot="none")
        ychartemp[3] <- chartemp$p
        ychartemp[4] <- chartemp$cv2
        ychartemp[5] <- fit[2]
      } else {
        ychartemp[6] = 0
      }
    } else {
      # Identify model
      # chartemp <- idclass(ytemp,type="PK",a.in=w.in,outplot="none")
      chartemp <- idclass(ytemp,type="PKa",a.in=w.in,outplot="none")
      ychartemp[3] <- chartemp$p
      ychartemp[4] <- chartemp$cv2
      ychartemp[5] <- which(chartemp$summary==1)
    }
  }
  
  # Fit model and produce forecasts
  if (ychartemp[6]==1){
    # Check if model.fit is given or optimise models
    if (!is.null(model.fit)){
      # Use model.fit
      if (fit[2] == 1){
        # Croston
        ftemp <- crost(ytemp,h=1,w=fit[3:4],init=fit[5:6],init.opt=FALSE)
      } else if (fit[2] == 2){
        # SBA
        ftemp <- crost(ytemp,h=1,w=fit[3:4],init=fit[5:6],init.opt=FALSE,type="sba")
      } else {
        # SES
        ftemp <- sexsm(ytemp,h=1,w=fit[3],init=fit[5],init.opt=FALSE)
      }
      p <- fit
    } else {
      # Fit new models
      if (ychartemp[5]==1){
        # Croston
        ftemp <- crost(ytemp,h=1,w=w,init.opt=init.opt)
        p <- c(minimumAL+i-1,1,ftemp$weights,ftemp$initial)
      } else if (ychartemp[5]==2){
        # SBA
        ftemp <- crost(ytemp,h=1,w=w,init.opt=init.opt,type="sba")
        p <- c(minimumAL+i-1,2,ftemp$weights,ftemp$initial)
      } else {
        ftemp <- sexsm(ytemp,h=1)
        p <- c(minimumAL+i-1,3,ftemp$alpha,NA,ftemp$initial,NA)
      }
    }
    # Disaggregate forecasts
    ftemp.out <- array(NA,c(1,h))
    ftemp.out[] <- ftemp$frc.out/ychartemp[1]
    fin <- as.vector(matrix(rep(ftemp$frc.in/ychartemp[1],ychartemp[1]),nrow=ychartemp[1],byrow=TRUE))
    ftemp.in <- array(NA,c(1,n))
    ftemp.in[yaggr$idx[[i]][1]:n] <- fin
    # f.out[i,] <- ftemp$frc.out/ychar[1,i]
    # f.in[i,(yaggr$idx[[i]][1]:n)] <- fin
  } else {
    # No model - not enough sample
    ftemp.out <- array(NA,c(1,h))
    ftemp.in <- array(NA,c(1,n))
    p <- c(minimumAL+i-1,rep(NA,5))
  }
  p <- matrix(p,ncol=6)

  return(cbind(ychartemp,p,ftemp.out,ftemp.in))

}