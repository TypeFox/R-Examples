#' Polar plot of temporal data
#' 
#' Representation in polar axis of the distribution of dates in the provided 
#' data set.
#' 
#' This function returns a plot representing the temporal distribution of 
#' records in the data set. This is done by representing dates in a radial axis,
#' with the distance from the center being the amount of records for that 
#' particular date. This function allows several arguments indicating different 
#' representation types. See the \code{arguments} section for an enumeration of 
#' them.
#' 
#' @import sqldf
#' @import plotrix
#' @importFrom grDevices col2rgb
#' @param indf input data frame containing biodiversity data set
#' @param timescale Temporal scale of the graph, or how are dates aggregated. 
#'   Accepted values are: d (daily, each feature in the plot represents a day), 
#'   w (weekly, each feature in the plot represents a week) and m (monthly, each
#'   feature in the plot represents a month). Default is d (daily).
#' @param title Title for the graph. Default is "Temporal coverage".
#' @param color color of the graph plot. Default is "red".
#' @param plottype Type of feature. Accepted values are: r (lines), p (polygon) 
#'   and s (symbols). Default is p (polygon).
#' @param avg If TRUE plots a graph of the average records rather than total 
#'   numbers. Default is FALSE.
#' @references Otegui, J., Arino, A. H., Encinas, M. A., & Pando, F. (2013). 
#'   Assessing the Primary Data Hosted by the Spanish Node of the Global 
#'   Biodiversity Information Facility (GBIF). PLoS ONE, 8(1), e55144. 
#'   doi:10.1371/journal.pone.0055144
#' @examples \dontrun{
#' tempolar(inat)
#' }
#' @family Temporal visualizations
#' @export
tempolar <- function(indf=NA, timescale=NA, title=NA, color=NA, plottype=NA,avg=FALSE){
  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)), 
               error = function(e) FALSE)
    })
  }
  if (!is.na(title)) {
    title2 <- title
  } else {
    title2 <- "Temporal coverage"
  }
  if (!is.na(color) & areColors(color)) {
    color2 <- color
  } else {
    color2 <- "red"
  }
  if (!is.na(plottype)) {
    plottype2 <- plottype
  } else {
    plottype2 <- "p"
  }
  if (!is.na(timescale)) {
    timescale2 <- timescale
  } else {
    timescale2 <- "d"
  }
  
  names(indf)=gsub("\\.","_",names(indf))
  if("Date_collected" %in% colnames(indf)){
    if(length(which(!is.na(indf$Date_collected)))==0){
      stop("Date_collected has no data")
    }
    dayofYear = as.numeric(strftime(as.Date(indf$Date_collected,na.rm=T), format = "%j"))
    weekofYear = as.numeric(strftime(as.Date(indf$Date_collected,na.rm=T), format = "%U"))
    monthofYear = as.numeric(strftime(as.Date(indf$Date_collected,na.rm=T), format = "%m"))
    Year_ = as.numeric(strftime(as.Date(indf$Date_collected,na.rm=T), format = "%Y"))
    
  } else {
    stop("Date_collected not found in data. Please use format_bdvis() to fix the problem")
  }
  indf = cbind(indf,dayofYear,weekofYear,monthofYear,Year_)
  if(timescale2=="d"){
    daytab=sqldf("select dayofYear, count(*) as dct from indf group by dayofYear")
    if(avg==F){
      if(is.na(daytab[1,1])){daytab=daytab[2:dim(daytab)[1],]}
      radial.plot(daytab$dct,
                  ((((daytab$dayofYear-1)*360)/366)*(3.14/180)),
                  line.col=color2, labels=month.abb,
                  clockwise=T, start=1.62,
                  radial.lim = c(0,max(daytab$dct)),
                  main=title2,boxed.radial=FALSE,
                  show.grid.labels=3,rp.type=plottype2)
    } else {
      alldays=sqldf("select dayofYear, Year_, count(*) as ct from indf group by dayofYear,monthofYear,Year_")
      daymean=sqldf("select dayofyear,avg(ct) as avgct,stdev(ct) as sdct from alldays group by dayofyear")
      
      if(is.na(daymean[1,1])){daymean=daymean[2:dim(daymean)[1],]}
      radial.plot(daymean$avgct,
                  ((((daymean$dayofYear-1)*360)/366)*(3.14/180)),
                  line.col=color2, labels=month.abb,
                  clockwise=T, start=1.62,
                  radial.lim = c(0,max(daymean$avgct)),
                  main=title2,boxed.radial=FALSE,
                  show.grid.labels=3,rp.type=plottype2)
    }
  }
  if(timescale2=="w"){
    weektab=sqldf("select weekofYear, count(*) as wct from indf group by weekofYear")
    if(avg==F){
      if(is.na(weektab[1,1])){weektab=weektab[2:dim(weektab)[1],]}
      if(dim(weektab)[1]==54){
        weektab[1,2]=weektab[1,2]+weektab[54,2]
        weektab=weektab[1:53,]
      }
      radial.plot(weektab$wct,
                  ((((weektab$weekofYear-1)*360)/53)*(3.14/180)),
                  line.col=color2,start=1.62, labels=month.abb,
                  radial.lim = c(0,max(weektab$wct)),
                  clockwise=TRUE,main=title2,boxed.radial=FALSE,
                  show.grid.labels=3,rp.type=plottype2,lwd=4)
    } else {
      allweeks=sqldf("select weekofYear, Year_, count(*) as ct from indf group by weekofYear,Year_")
      weekmean=sqldf("select weekofyear,avg(ct) as avgct,stdev(ct) as sdct from allweeks group by weekofyear")
      if(is.na(weekmean[1,1])){weekmean=weekmean[2:dim(weekmean)[1],]}
      if(dim(weekmean)[1]==54){
        weekmean[1,2]=weekmean[1,2]+weekmean[54,2]
        weekmean=weekmean[1:53,]
      }
      radial.plot(weekmean$avgct,
                  ((((weekmean$weekofYear-1)*360)/53)*(3.14/180)),
                  line.col=color2,start=1.62, labels=month.abb,
                  radial.lim = c(0,max(weekmean$avgct)),
                  clockwise=TRUE,main=title2,boxed.radial=FALSE,
                  show.grid.labels=3,rp.type=plottype2,lwd=4)
    }
  }
  if(timescale2=="m"){
    monthtab=sqldf("select monthofYear, count(*) as mct from indf group by monthofYear")
    if(avg==F){
      if(is.na(monthtab[1,1])){monthtab=monthtab[2:dim(monthtab)[1],]}
      radial.plot(monthtab$mct,
                  ((((monthtab$monthofYear-1)*360)/12)*(3.14/180)),
                  line.col=color2,start=1.62, labels=month.abb,
                  radial.lim = c(0,max(monthtab$mct)),
                  clockwise=TRUE,main=title2,boxed.radial=FALSE,
                  show.grid.labels=3,rp.type=plottype2,lwd=4)  
    } else {
      allmonths=sqldf("select monthofYear, Year_, count(*) as ct from indf group by monthofYear,Year_")
      monthmean=sqldf("select monthofyear,avg(ct) as avgct,stdev(ct) as sdct from allmonths group by monthofyear")
      if(is.na(monthmean[1,1])){monthmean=monthmean[2:dim(monthmean)[1],]}
      radial.plot(monthmean$avgct,
                  ((((monthmean$monthofYear-1)*360)/12)*(3.14/180)),
                  line.col=color2,start=1.62, labels=month.abb,
                  radial.lim = c(0,max(monthmean$avgct)),
                  clockwise=TRUE,main=title2,boxed.radial=FALSE,
                  show.grid.labels=3,rp.type=plottype2,lwd=4)  
    }
  }
}