#' Spline-based window analysis
#'
#' Defines window boundaries for analyzing genetic data and evaluates the specified windows
#'
#' @param Y A vector of estimates of some parameter, such as Fst, computed at individual markers. One chromosome should be analyzed at a time.
#' @param map A vector of positions for each marker with a corresponding parameter estimate.
#' @param smoothness The level of resolution (in base pairs) for computing the spline and its derivatives
#' @param s2 The variance of parameter estimates, to be used for computing Wstats. Default is to compute this automatically, but it may be manually specified, e.g. so that the value across chromosomes may be utilized.
#' @param mean The mean of parameter estimates, to be used for computing Wstats. Default is to comput this automatically, but it may be manually specified, e.g. so that the value across chromosomes may be utilized.
#' @param plotRaw Whether or not to produce a plot of raw data, with the fitted spline
#' @param plotWindows Whether or not to include a plot of Wstat values over the computed windows
#' @param method The method for controlling amount of smoothing: 1, 2, 3, or 4. See documentation of smooth.Pspline for description. Usual choices are either 3 for generalized cross validation or 4 for ordinary cross validation.
#' 
#' @return rawSpline The fitted spline object
#' @return breaks The spline-suggested window breaks
#' @return windowData A table of mean parameter estimates and Wstats computed over spline-suggested windows
#'
#' @examples
#' 
#' data(chr6)
#' sub6 <- chr6[55000:63000,]
#' chr6Spline <- splineAnalyze(Y=sub6$Fst,map=sub6$Position,smoothness=100,
#' plotRaw=TRUE,plotWindows=TRUE,method=4)
#' 
#'\dontrun{
#' chr6Spline <- splineAnalyze(Y=chr6$Fst,map=chr6$Position,smoothness=100,
#' plotRaw=TRUE,plotWindows=TRUE,method=4)
#'}
#'
#' @import pspline
#' @export



splineAnalyze <- function(Y,map,smoothness=100,s2=NA,mean=NA,plotRaw=FALSE,plotWindows=FALSE, method=3){
  ### Function to find roots
  roots <- function(data){
   data1 <- c(NA,data[1:{length(data)-1}])
   data2 <- data
   posneg <- which(data1>0 & data2<0) - 0.5
   negpos <- which(data1<0 & data2>0) -0.5
   zero <- which(data == 0)
   roots <- sort(c(posneg,negpos,zero))
   return(roots)
  }
  #set up data sets
  rawData <- data.frame(Pos=map,Y=Y)
  data <- rawData[which(is.na(rawData$Y)==FALSE),]
  #compute spline and derivatives
  pspline <- smooth.Pspline(data[,1],data[,2],norder=2,method=method)
  predict <- predict(pspline,seq(0,max(pspline$x),by=smoothness))
  psplinederiv <- predict(pspline,seq(0,max(pspline$x),by=smoothness),nderiv=2)
  psplineInflection <- roots(psplinederiv)*smoothness
  # Print number of windows
  print(paste("Total number of windows = ", length(psplineInflection) + 1))
  #create window table
  print(" ---- Computing window statistics ----")
  if(is.na(s2)) s2 <- var(Y)
  if(is.na(mean)) mean <- mean(Y,na.rm=TRUE)
  cat(1, "of", length(psplineInflection)+1, "\r") # print progress
  Distinct <- data.frame(WindowStart=rep(NA,length(psplineInflection)+1),WindowStop=NA,SNPcount=NA,MeanY=NA, Wstat=NA)
  Distinct$WindowStart[1] <- min(pspline$x)
  Distinct$WindowStop[1] <- psplineInflection[1]
  Distinct$SNPcount[1] <- length(which(data[,1] <= psplineInflection[1]))
  Distinct$MeanY[1] <- mean(data[which(data[,1] <= psplineInflection[1]),2],na.rm=TRUE)
  Distinct$Wstat[1] <- {mean(data[which(data[,1] <= psplineInflection[1]),2],na.rm=TRUE)-mean}/
                        sqrt(s2/length(which(data[,1] <= psplineInflection[1])))
  for(i in 2:length(psplineInflection)){
    cat(i, "of", length(psplineInflection)+1, "\r")
    Distinct$WindowStart[i] <- psplineInflection[i-1]
    Distinct$WindowStop[i] <- psplineInflection[i]
    Distinct$SNPcount[i] <- length(which(data[,1] >= psplineInflection[i-1] & data[,1]
                                   <= psplineInflection[i]))
    Distinct$MeanY[i] <- mean(data[which(data[,1] >= psplineInflection[i-1] & data[,1]
                           <= psplineInflection[i]),2],na.rm=TRUE)
    Distinct$Wstat[i] <- {mean(data[which(data[,1] >= psplineInflection[i-1] & data[,1]
                         <= psplineInflection[i]),2],na.rm=TRUE) - mean}/ sqrt(s2/length(which(data[,1]
                         >= psplineInflection[i-1] & data[,1] <= psplineInflection[i])))
  }
  ### Fill out final window
  i <- i+1
  cat(i, "of", length(psplineInflection)+1, "\r")
  Distinct$WindowStart[i] <- psplineInflection[i-1]
  Distinct$WindowStop[i] <- max(data$Pos)
  Distinct$SNPcount[i] <- length(which(data[,1] >= psplineInflection[i-1] ))
  Distinct$MeanY[i] <- mean(data[which(data[,1] >= psplineInflection[i-1]),2],na.rm=TRUE)
  Distinct$Wstat[i] <- {mean(data[which(data[,1] >= psplineInflection[i-1]),2],na.rm=TRUE) - mean}/
                        sqrt(s2/length(which(data[,1] >= psplineInflection[i-1])))

  print(" ---- done ---- ")
  #make plots if requested
  if(plotRaw==TRUE & plotWindows==TRUE) par(mfrow=c(2,1))
  if(plotRaw==TRUE){
    plot(data,xlab="Position (bp)",ylab="Raw values")
    lines(seq(0,max(pspline$x),by=smoothness),predict,col="red")
  }
  if(plotWindows==TRUE){
    plot((Distinct$WindowStop-Distinct$WindowStart)/2+Distinct$WindowStart,Distinct$Wstat,xlab="Position (bp)", ylab="Spline Wstat",pch=19)
  }
  return(list(rawSpline=pspline,breaks=psplineInflection,windowData=Distinct))
}
