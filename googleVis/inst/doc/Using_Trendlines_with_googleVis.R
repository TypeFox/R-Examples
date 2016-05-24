## ----setOptions, message=FALSE, echo=FALSE-------------------------------
library(googleVis)
library(knitr)
op <- options(gvis.plot.tag='chart')
read_demo('Trendlines', 'googleVis')

## ----LinearTrend, results='asis', tidy=FALSE-----------------------------
plot(
  gvisScatterChart(women, options=list(trendlines="0"))
)

## ----ExponentialTrend, results='asis', tidy=FALSE------------------------
plot(
  gvisScatterChart(women, options=list(
    trendlines="{0: { type: 'exponential',  
                     visibleInLegend: 'true', 
                     color: 'green',
                     lineWidth: 10,
                     opacity: 0.5}}",
    chartArea="{left:50,top:20,width:'50%',height:'75%'}"))
)

## ----ColumnChartWithTrendline, results='asis', tidy=FALSE----------------
dat <- data.frame(val1=c(1,3,4,5,6,8), 
                  val2=c(12,23,32,40,50,55))
plot(
  gvisColumnChart(dat,
                  options=list(trendlines="{0: {}}"))
)

## ----DifferentLabels, results='asis', tidy=FALSE-------------------------
dat$val3 <- c(5,6,10,12,15,20)
plot(
  gvisColumnChart(dat,
                  options=list(trendlines="{
                          0: {
                            labelInLegend: 'Trendline 1',
                            visibleInLegend: true,}, 
                          1:{
                            labelInLegend: 'Trendline 2',
                            visibleInLegend: true}
                          }",
                          chartArea="{left:50,top:20,
                                      width:'50%',height:'75%'}"
                  ))
)

