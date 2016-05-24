#'Draws a chronohorogram of records
#'
#'Draws a detailed temporal representation (also known as chronohorogram) of the
#'dates in the provided data set. For more information on the chronhorogram, 
#'please see the \code{References} section.
#'
#'@import sqldf
#'@import plotrix
#'@importFrom grDevices colorRampPalette
#'@param indf input data frame containing biodiversity data set
#'@param title title of the plot. Default is "Chronohorogram"
#'@param startyear starting year for the plot. Default is 1980
#'@param endyear end year for the plot. Default is current year
#'@param colors Pair of colors to build color gradient, in the form of a
#'  character vector. Default is blue (less) - red (more) gradient
#'  \code{c("red", "blue")}
#'@references Arino, A. H., & Otegui, J. (2008). Sampling biodiversity sampling.
#'  In Proceedings of TDWG (pp. 77-78). Retrieved from 
#'  http://www.tdwg.org/fileadmin/2008conference/documents/Proceedings2008.pdf#page=77
#'  
#'  
#'@family Temporal visualizations
#'@export
#'@examples \dontrun{
#'chronohorogram(inat)
#'}
chronohorogram <- function (indf=NA,title="Chronohorogram",startyear=1980,endyear=NA,
                            colors=c("red", "blue")){
  if(!is.element("Date_collected",names(indf))){
    stop("Field Date_collected not found in dataset \n")
  }
  if (is.na(endyear)){
    endyear=as.integer(format(Sys.Date(), "%Y"))
  }
  dat1=sqldf("select Date_collected, count(*) as ct from indf group by Date_collected")
  if(is.na(dat1$Date_collected[1])){dat1=dat1[2:dim(dat1)[1],]}
  if(as.character(dat1$Date_collected[1])==""){dat1=dat1[2:dim(dat1)[1],]}
  d=as.numeric(strftime(as.Date(dat1$Date_collected,na.rm=T), format = "%j"))
  y=as.numeric(strftime(as.Date(dat1$Date_collected,na.rm=T), format = "%Y"))
  rind=which(y<=endyear & y>=startyear)
  d=d[rind]
  y=y[rind]
  crp=colorRampPalette(colors)( round(log10(max(dat1$ct)+2)) +1 ) 
  radial.plot(y,(d/366)*360,rp.type="s", start=1.62, labels=month.abb,
              clockwise=TRUE, point.col=crp[(round(log10(dat1[,2])+2))-1], 
              point.symbols=20, grid.bg="black",radial.lim=c(startyear,endyear),
              show.radial.grid=T,show.grid.labels=F,show.grid=T,
              main=title)
}