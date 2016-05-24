NULL
#' 
#' Plots curves for comparison between measures and simulations (hourly interpolations)
#'
#' @title Plots hourly curves
#'
#' @author  Emanuele Eccel, Emanuele Cordano \email{emanuele.eccel@@iasma.it}
#' 
#' @param meas  measured hourly values file (table), where the first column is the series' ID
#' @param sim  simulated hourly list (output of Th_int_series)
#' @param series  names of the seris to plot. If \code{NULL}, plots all the series
#' @param chart.start  start date for the plotting. Format example: "1Jan2000"
#' @param chart.end  end date for the plotting. Format example: "1Jan2000"
#' @param date.format input date format for measurements (formats for function \code{chron}). Default is "ymd"
#' @param missing_code  code (either real or character) for missing values in measurements. Default is \code{NA}
#' @param wait  lag time (seconds) between plot appearance on the screen (default is 1 second)
#' @param plot.leg  logical: if \code{TRUE} (default) legends are plotted
#' @param leg.pos  position of legends (only if \code{plot.leg} = \code{TRUE}). Default is "bottomright". Values for \code{par} can be passed to function

#' @export 

#' @return A plot with two curves: measured values and hourly interpolations

#' @references 
#' Eccel, E., 2010: What we can ask to hourly temperature recording. Part II: hourly interpolation of temperatures for climatology and modelling. Italian Journal of Agrometeorology XV(2):45-50 
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag45.pdf},\url{www.agrometeorologia.it}
#' 
#' 
#' See also: Eccel, E., 2010: What we can ask to hourly temperature recording. Part I: statistical vs. meteorological meaning of minimum temperature. Italian Journal of Agrometeorology XV(2):41-43.
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag41.pdf}
#' 
#' @note
#' If daily minimum and maximum values are the absolute ones for each day, the interpolated curve will be generally either higher (maximum) or lower (minimum) than hourly measurements, which are mean hourly values (hence, slightly lower and higher than the absolute daily ones, respectively). This may explain an apparent mismatch between the two curves.
#'

#' @examples
#' data(Trentino_hourly_T)
#' stations <- c("T0001","T0010","T0129")
#' plot_meas_sim(meas = h_d_t, sim = Th_int_list, series=stations,
#'               missing_code = -999.9, chart.start = "1Feb2004", 
#'              chart.end = "29Feb2004", leg.pos = "top")


########################################################
# PLOTS SERIES OF MEASURED AND SIMULATED TEMPERATURES
########################################################


plot_meas_sim<-function(meas, sim, series=NULL, chart.start, chart.end, date.format="ymd", missing_code=NA, wait=1, plot.leg=TRUE, leg.pos="bottomright")

{

if(is.null(series)) series<-unique(meas[,1]) 

for(sta in series)
{

# selects and filters measures
T_meas_hourly<-meas[meas[,1]==sta,] 
data.date<-as.date(as.character(T_meas_hourly[,2]), order=date.format)
data<-date.mdy(data.date)
T_meas_hourly[,4][T_meas_hourly[,4]==missing_code]<-NA
names(T_meas_hourly)<-c("ID", "date", "hour","T")
if(min(data.date) > as.date(chart.start) | max(data.date) < as.date(chart.end))
  print(paste(sta, "Warning: required period of instrumental recording exceeds data availability"), quote=FALSE)
T_meas_x.graph<-T_meas_hourly[data.date>= as.date(chart.start) & data.date<= as.date(chart.end),]

# selects and filters simulations
data.date<-as.date(paste(sim$Date$year,sim$Date$month, sim$Date$day, sep="/"), order=date.format)
if(min(data.date) > as.date(chart.start) | max(data.date) < as.date(chart.end))
  print(paste(sta, "Warning: required simulation period exceeds data availability"), quote=FALSE)
T_sim_x.graph<-sim[[sta]][data.date>= as.date(chart.start) & data.date<= as.date(chart.end)]

if(sum(!is.na(T_sim_x.graph))>0 & sum(!is.na(T_meas_x.graph$T))>0 )  # there are data
{
ylim<-c(min(c(T_meas_x.graph$T,T_sim_x.graph), na.rm=TRUE), max(c(T_meas_x.graph$T,T_sim_x.graph), na.rm=TRUE))
plot(T_meas_x.graph$T, type="l", ylim=ylim, main=paste(sta,"    ",chart.start,"-", chart.end), xaxp=c(0,nrow(T_meas_x.graph),as.date(chart.end)-as.date(chart.start)+1), ylab="T (deg C)", xlab="hours from start")
par(new=TRUE)
plot(T_sim_x.graph, type="l", ylim=ylim, xaxp=c(0,nrow(T_meas_x.graph),as.date(chart.end)-as.date(chart.start)+1), ylab="", xlab="", lty="dashed", col=2)
if(plot.leg==TRUE) legend(leg.pos, cex=0.6, c("meas.", "simul."), lty=c("solid", "dashed"), col=c(1,2), bty="o", bg=NULL)
Sys.sleep(wait)
}

}  # sta

}