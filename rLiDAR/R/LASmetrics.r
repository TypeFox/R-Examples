#'LiDAR metrics computation
#'
#'@description Compute statistically based metrics that describe the LiDAR dataset
#'
#'@usage LASmetrics(LASfile, minht, above)
#'
#'@param LASfile A LAS standard LiDAR data file
#'@param minht Use only returns above specified height break, e.g. 1.30 m. Default is 1.37 m.
#'@param above Compute covers metrics using specified height break, e.g. 2.5 m. Default is 2 m.
#'@return Returns A matrix with the LiDAR-derived vegetation height and canopy cover metrics (see \emph{cloudmetrics}, in McGaughey, 2014)
#'@author Carlos Alberto Silva
#'@seealso McGaughey, R. 2014. FUSION/LDV: Software for lidar data analysis and visualization. Version 3.41. Seattle, WA: U.S. Department of Agriculture, Forest Service, Pacific Northwest Research Station.
#'@examples
#'
#'#=======================================================================#
#'# Example 01: Computing LiDAR metrics for a single LAS file
#'#=======================================================================#
#'# Import the LAS data file
#'LASfile <- system.file("extdata", "LASexample1.las", package="rLiDAR")
#'
#'# Set the minht and above parameters
#'minht<-1.37  # meters or feet
#'above<-2.00  # meters or feet
#'
#'# LiDAR metrics computation
#'LiDARmetrics<-LASmetrics(LASfile, minht, above)
#'
#'#==========================================================================#
#'# Example 02: Computing LiDAR metrics for multiple LAS files within a folder
#'#==========================================================================#
#'# Set folder where LAS source files reside
#'folder=dirname(LASfile)
#'
#'# Get list of LAS files residing in the folder
#'LASlist <- list.files(folder, pattern="*.las", full.names=TRUE)
#'
#'# Set the "minht" and "above" parameters
#'minht<-1.37  # meters or feet
#'above<-2.00  # meters or feet
#'
#'# Creat an empty dataframe in whic to store the LiDAR metrics
#'getMetrics<-data.frame()
#'
#'# Set a loop to compute the LiDAR metrics
#'for ( i in LASlist) {
#'  getMetrics<-rbind(getMetrics, LASmetrics(i, minht, above))}
#'
#'# Table of the LiDAR metrics
#'LiDARmetrics<-cbind(Files=c(basename(LASlist)), getMetrics)
#'head(LiDARmetrics)
#'
#'@importFrom bitops bitAnd bitShiftR
#'@export
LASmetrics<-function(LASfile,minht=1.37,above=2) {

  if (class(minht)!="numeric") {stop("The minht parameter is invalid. It is not a numeric input")}
  if (class(above)!="numeric") {stop("The above parameter is invalid. It is not a numeric input")}
  
  LASfile<-readLAS(LASfile, short=T)
  
  MaxZ<-max(LASfile[,"Z"])
  
  if (minht >= MaxZ) {stop(paste0("The minht parameter is invalid. It must to be less than ",MaxZ))}
  
  allreturn<-nrow(LASfile)
  allreturn_minht<-subset(LASfile,LASfile[,"Z"] > minht)
  firstReturn<-subset(LASfile, LASfile[,"ReturnNumber"] == 1)
  allreturnAbove<-subset(LASfile, LASfile[,"Z"] > above)
  firstReturnAbove<-subset(firstReturn, firstReturn[,"Z"] > above)
  firstabovemean<-subset(firstReturn, firstReturn[,"Z"] > mean(allreturnAbove[,"Z"]))
  firstabovemode<-subset(firstReturn, firstReturn[,"Z"] > as.numeric(names(table(allreturnAbove[,"Z"]))[which.max(table(allreturnAbove[,"Z"]))]))
  allabovemean<-subset(LASfile, LASfile[,"Z"] > mean(allreturnAbove[,"Z"]))
  hmode<-(names(table(allreturnAbove[,"Z"]))[which.max(table(allreturnAbove[,"Z"]))])
  allabovemode<-subset(LASfile, LASfile[,"Z"] > hmode)
  
  # from moments package: Lukasz Komsta et al.  (2015) -----------------# 
  "skewness" <-
    function (x, na.rm = FALSE) 
    {
      if (is.matrix(x)) 
        apply(x, 2, skewness, na.rm = na.rm)
      else if (is.vector(x)) {
        if (na.rm) x <- x[!is.na(x)] 
        n <- length(x)
        (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
      }
      else if (is.data.frame(x)) 
        sapply(x, skewness, na.rm = na.rm)
      else skewness(as.vector(x), na.rm = na.rm)
    }
  
  "kurtosis" <-
    function (x, na.rm = FALSE) 
    {
      if (is.matrix(x)) 
        apply(x, 2, kurtosis, na.rm = na.rm)
      else if (is.vector(x)) {
        if (na.rm) x <- x[!is.na(x)] 
        n <- length(x)
        n*sum( (x-mean(x))^4 )/(sum( (x-mean(x))^2 )^2)
      }
      else if (is.data.frame(x)) 
        sapply(x, kurtosis, na.rm = na.rm)
      else kurtosis(as.vector(x), na.rm = na.rm)
    }
  #------------------------------------------------------------------#
  metrics<-data.frame(

  Total.all.return.count=nrow(LASfile),
    Total.first.return.count=nrow(firstReturn),
    Total.all.return.count.aboveXX=nrow(allreturn_minht),
    
    Return.1.count=nrow(subset(allreturn_minht,allreturn_minht[,"ReturnNumber"]==1)),
    Return.2.count=nrow(subset(allreturn_minht, allreturn_minht[,"ReturnNumber"]==2)),
    Return.3.count=nrow(subset(allreturn_minht, allreturn_minht[,"ReturnNumber"]==3)),
    Return.4.count=nrow(subset(allreturn_minht, allreturn_minht[,"ReturnNumber"]==4)),
    Return.5.count=nrow(subset(allreturn_minht, allreturn_minht[,"ReturnNumber"]==5)),
    Return.6.count=nrow(subset(allreturn_minht, allreturn_minht[,"ReturnNumber"]==6)),
    Return.7.count=nrow(subset(allreturn_minht, allreturn_minht[,"ReturnNumber"]==7)),
    Return.8.count=nrow(subset(allreturn_minht, allreturn_minht[,"ReturnNumber"]==8)),
    Return.9.count=nrow(subset(allreturn_minht, allreturn_minht[,"ReturnNumber"]==9)),
    
    hmin=round(min(allreturn_minht[,"Z"]), digits=2),
    hmax=round(max(allreturn_minht[,"Z"]), digits=2),
    hmean=round(mean(allreturn_minht[,"Z"]),digits=2),
    hmode = round(as.numeric(names(table(allreturn_minht[,"Z"]))[which.max(table(allreturn_minht[,"Z"]))]), digits=2),
    
    hmedian=round(median(allreturn_minht[,"Z"]),digits=2),
    hsd=round(sd(allreturn_minht[,"Z"]),digits=2),
    hvar=round(var(allreturn_minht[,"Z"]),digits=2),
    hcv=round((sd(allreturn_minht[,"Z"])/mean(allreturn_minht[,"Z"]))*100,digits=2),
    hkurtosis=round(kurtosis(allreturn_minht[,"Z"]),digits=2),
    hskewness=round(skewness(allreturn_minht[,"Z"]),digits=2),
    hP1=round(quantile(allreturn_minht[,"Z"],0.01),digits=2),
    hP5=round(quantile(allreturn_minht[,"Z"],0.05),digits=2),
    hP10=round(quantile(allreturn_minht[,"Z"],0.1),digits=2),
    hP15=round(quantile(allreturn_minht[,"Z"],0.15),digits=2),
    hP20=round(quantile(allreturn_minht[,"Z"],0.20),digits=2),
    hP25=round(quantile(allreturn_minht[,"Z"],0.25),digits=2),
    hP30=round(quantile(allreturn_minht[,"Z"],0.30),digits=2),
    hP35=round(quantile(allreturn_minht[,"Z"],0.35),digits=2),
    hP40=round(quantile(allreturn_minht[,"Z"],0.40),digits=2),
    hP45=round(quantile(allreturn_minht[,"Z"],0.45),digits=2),
    hP50=round(quantile(allreturn_minht[,"Z"],0.50),digits=2),
    hP55=round(quantile(allreturn_minht[,"Z"],0.55),digits=2),
    hP60=round(quantile(allreturn_minht[,"Z"],0.60),digits=2),
    hP65=round(quantile(allreturn_minht[,"Z"],0.65),digits=2),
    hP70=round(quantile(allreturn_minht[,"Z"],0.70),digits=2),
    hP75=round(quantile(allreturn_minht[,"Z"],0.75),digits=2),
    hP80=round(quantile(allreturn_minht[,"Z"],0.85),digits=2),
    hP90=round(quantile(allreturn_minht[,"Z"],0.90),digits=2),
    hP95=round(quantile(allreturn_minht[,"Z"],0.95),digits=2),
    hP99=round(quantile(allreturn_minht[,"Z"],0.99),digits=2),
       
    Canopy.relief.ratio= ((mean(allreturn_minht[,"Z"])-min(allreturn_minht[,"Z"]))/(max(allreturn_minht[,"Z"])-min(allreturn_minht[,"Z"]))),
    
    Imin=round(min(allreturn_minht[,"Intensity"],digits=2)),
    Imax=round(max(allreturn_minht[,"Intensity"]),digits=2),
    Imean=round(mean(allreturn_minht[,"Intensity"]),digits=2),
    Imode = round(as.numeric(names(table(allreturn_minht[,"Intensity"]))[which.max(table(allreturn_minht[,"Intensity"]))]),digits=2),
    Imedian=round(median(allreturn_minht[,"Intensity"]),digits=2),
    Isd=round(sd(allreturn_minht[,"Intensity"]),digits=2), 
    Ivar=round(var(allreturn_minht[,"Intensity"]),digits=2),
    Icv=round((sd(allreturn_minht[,"Intensity"])/mean(allreturn_minht[,"Intensity"]))*100,digits=2),
    Ikurtosis=round(kurtosis(allreturn_minht[,"Intensity"]),digits=2),
    Iskewness=round(skewness(allreturn_minht[,"Intensity"]),digits=2),
    IP1=round(quantile(allreturn_minht[,"Intensity"],0.01),digits=2),
    IP5=round(quantile(allreturn_minht[,"Intensity"],0.05),digits=2),
    IP10=round(quantile(allreturn_minht[,"Intensity"],0.1),digits=2),
    IP15=round(quantile(allreturn_minht[,"Intensity"],0.15),digits=2),
    IP20=round(quantile(allreturn_minht[,"Intensity"],0.20),digits=2),
    IP25=round(quantile(allreturn_minht[,"Intensity"],0.25),digits=2),
    IP30=round(quantile(allreturn_minht[,"Intensity"],0.30),digits=2),
    IP35=round(quantile(allreturn_minht[,"Intensity"],0.35),digits=2),
    IP40=round(quantile(allreturn_minht[,"Intensity"],0.40),digits=2),
    IP45=round(quantile(allreturn_minht[,"Intensity"],0.45),digits=2),
    IP50=round(quantile(allreturn_minht[,"Intensity"],0.50),digits=2),
    IP55=round(quantile(allreturn_minht[,"Intensity"],0.55),digits=2),
    IP60=round(quantile(allreturn_minht[,"Intensity"],0.60),digits=2),
    IP65=round(quantile(allreturn_minht[,"Intensity"],0.65),digits=2),
    IP70=round(quantile(allreturn_minht[,"Intensity"],0.70),digits=2),
    IP75=round(quantile(allreturn_minht[,"Intensity"],0.75),digits=2),
    IP80=round(quantile(allreturn_minht[,"Intensity"],0.85),digits=2),
    IP90=round(quantile(allreturn_minht[,"Intensity"],0.90),digits=2),
    IP95=round(quantile(allreturn_minht[,"Intensity"],0.95),digits=2),
    IP99=round(quantile(allreturn_minht[,"Intensity"],0.99),digits=2),
    
    Pentage.first.returns.Above.XX=(nrow(firstReturnAbove))/(nrow(firstReturn))*100,
    Percentage.all.returns.above.XX=(nrow(allreturnAbove)/allreturn)*100,
    All.returns.above.XX.Total.first.returns.100=(nrow(allreturnAbove)/nrow(firstReturn))*100,
    First.returns.above.XX=nrow(firstReturnAbove),
    All.returns.above.XX=nrow(allreturnAbove),
    
    Percentage.first.returns.above.mean = (nrow(firstabovemean)/nrow(firstReturn))*100,
    Percentage.first.returns.above.mode = (nrow(firstabovemode)/nrow(firstReturn))*100,
    Percentage.all.returns.above.mean = (nrow(allabovemean)/allreturn)*100,
    Percentage.all.returns.above.mode =  (nrow(allabovemode)/allreturn)*100,
    
    All.returns.above.mean.Total.first.returns.100 = (nrow(allabovemean)/nrow(firstReturn))*100,
    All.returns.above.mode.Total.first.returns.100 = (nrow(allabovemode)/nrow(firstReturn))*100,
    First.returns.above.mean= nrow(firstabovemean),
    First.returns.above.mode= nrow(firstabovemode),
    All.returns.above.mean= nrow(allabovemean),
    All.returns.above.mode= nrow(allabovemode))
  
  
  
  colnames(metrics)<-c("Total all return count","Total first return count",paste("Total all return count", above),paste("Return 1 count", above),paste("Return 2 count", above),paste("Return 3 count", above),paste("Return 4 count", above),
                       paste("Return 5 count", above),paste("Return 6 count", above),paste("Return 7 count", above),paste("Return 8 count", above),paste("Return 9 count", above),"hmin",
                       "hmax","hmean",  "hmode","hmedian","hsd","hvar","hcv","hkurtosis","hskewness","hP1",
                       "hP5","hP10","hP15","hP20","hP25","hP30","hP35","hP40","hP45","hP50","hP55","hP60","hP65","hP70","hP75",
                       "hP80","hP90","hP95","hP99","Canopy.relief.ratio","Imin","Imax","Imean","Imode","Imedian","Isd","Ivar","Icv",
                       "Ikurtosis","Iskewness","IP1","IP5","IP10","IP15","IP20","IP25","IP30","IP35","IP40","IP45","IP50","IP55","IP60",
                       "IP65","IP70","IP75","IP80","IP90","IP95","IP99",paste("Pentage first returns Above",above),paste("Percentage all returns above",above),
                       paste("(All returns above",above,"/ Total first returns)*100"),paste("First returns above",above),paste("All returns above",above),"Percentage first returns above mean",
                       "Percentage first returns above mode","Percentage.all.returns.above.mean","Percentage all returns above mode","(All returns above mean / Total first returns)*100",
                       "(All returns above mode / Total first returns)* 100","First returns above mean","First returns above mode","All returns above mean","All returns above mode")
  rownames(metrics)<-NULL
  return(data.frame(metrics))
}

