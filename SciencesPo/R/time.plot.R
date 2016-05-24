if (getRversion() >= "2.15.1") globalVariables(c("value", "z", "time", "type", "lag", "acf", "pacf"))

#' @title Make the ggplot2 version of TS plots
#'
#' @description The function produces TS plots using ggplot.
#'
#' @param ts The TS object.
#' @param ylab The \code{y-axis} title.
#' @param ylim The \code{y-axis} limits.
#' @param ci The desired confidence interval.
#' @param \dots Ignored parameters passed to ggplot.
#'
#' @examples
#'  ts.sim <- stats::arima.sim(n = 100, list(ma=0.8), innov=rnorm(100))
#'  timeplot(ts.sim)
#' @export
#'
#' @import ggplot2
`timeplot` <-function(ts, ylab = '', ylim=c(-1,1), ci=.95, ...){
    ts.df <- data.frame(time=c(1:length(ts)),value=ts)
    timeSeriesPlot <- ggplot(ts.df,aes(x=time,y=value)) +
    geom_line()
    ts.acf <-stats::acf(ts, plot=FALSE)
    ts.pacf <-stats::pacf(ts, plot=FALSE)
    clim0 <- stats::qnorm((1 + ci)/2)/sqrt(ts.acf$n.used)
    clim <- c(-clim0,clim0)

    hline.data <- data.frame(z=c(0,clim),
                             type=c("base","ci","ci"))

acfPlot <- ggplot(data.frame(lag=ts.acf$lag,acf=ts.acf$acf)) +
 geom_hline(aes(yintercept=z,colour=type,linetype=type),hline.data) +
 geom_linerange(aes(x=lag,ymin=0,ymax=acf)) +
  scale_colour_manual(values = c("black","black")) +
  scale_linetype_manual(values =c("solid","dashed")) +
  ggtitle("Autocorrelations")

pacfPlot <- ggplot(data.frame(lag=ts.pacf$lag, pacf=ts.pacf$acf)) +
geom_hline(aes(yintercept=z,colour=type,linetype=type),hline.data) +
      geom_linerange(aes(x=lag,ymin=0,ymax=pacf)) +
      scale_colour_manual(values = c("black","black")) +
      scale_linetype_manual(values =c("solid","dashed")) +
      ggtitle("Partial Autocorrelations")
gridExtra::grid.arrange(timeSeriesPlot, gridExtra::arrangeGrob(acfPlot, pacfPlot, ncol=2),ncol=1)
}
