 require("tfplot")
 
# See also
# - blog post :  http://www.packtpub.com/article/creating-time-series-charts-r
# - chart.TimeSeries in xts?

z <- ts(rnorm(30), start=c(1990,1), frequency=12)
tfplot(z, Xaxis="auto", xlab= "time, time, time", source="Source: RNG")

tfplot(z, Xaxis="auto", xlab= "time, time, time", source="Source: RNG",
  YaxisR=TRUE, Yaxis.lab.rot = "horizontal")
tfplot(z, Xaxis="auto", xlab= "time, time, time", source="Source: RNG",
  YaxisL=FALSE,YaxisR=TRUE, Yaxis.lab.rot = "horizontal")

tfplot(z, Xaxis="auto", xlab= "time, time, time", source="Source: RNG")
tfplot(z, Xaxis="auto", xlab= "time, time, time", source="Source: RNG")

z <- ts(rnorm(50), start=c(1990,1), frequency=12)
tfplot(z, Xaxis=NULL,   xlab= "time, time, time", source="Source: RNG")
tfplot(z, Xaxis="auto", xlab= "time, time, time", source="Source: RNG")
tfplot(z, Xaxis="auto", xlab= "time, time, time", source="Source: RNG", start=c(1991,1))

z <- ts(matrix(rnorm(80),40,2), start=c(1990,1), frequency=12)
tfplot(z, Xaxis=NULL,   xlab= "time, time, time")
tfplot(z, Xaxis="auto", xlab= "time, time, time")

z <- ts(matrix(rnorm(24),12,2), start=c(1990,1), frequency=12)
tfplot(z, Xaxis="auto", xlab= "time, time, time")

z <- ts(rnorm(40), start=c(1990,6), frequency=12)
tfplot(z, Xaxis="auto", xlab= "time, time, time")

z <- ts(rnorm(10), start=c(1990,1), frequency=12)
tfplot(z, Xaxis="auto", xlab= "time, time, time")


z <- ts(rnorm(30), start=c(1990,2), frequency=4)
tfplot(z, Xaxis="auto", xlab= "time, time, time")

z <- ts(rnorm(40), start=c(1990,1), frequency=4)
tfplot(z, Xaxis="auto", xlab= "time, time, time")


z <- ts(rnorm(10), start=c(1990,1), frequency=1)
tfplot(z, Xaxis="auto", xlab= "time, time, time")
tfplot(z, Xaxis=NULL,   xlab= "time, time, time")

z <- ts(rnorm(40), start=c(1990,1), frequency=1)
tfplot(z, Xaxis="auto", xlab= "time, time, time")

# require(zoo)
# z <- zooreg(rnorm(40), start=c(1990,1), frequency=12)
# tfplot(z, xlab= "time, time, time")
# tfplot(z, Xaxis=NULL,   xlab= "time, time, time", source="Source: RNG")
# tfplot(z, Xaxis="auto", xlab= "time, time, time", source="Source: RNG")
# 
# z <- zoo(rnorm(100), order.by=as.Date("2003-02-01")+1:100)
# tfplot(z, xlab= "time, time, time", source="Source: RNG")
# tfplot(z, Xaxis="auto", xlab= "time, time, time", source="Source: RNG")
# 
# z <- zoo(rnorm(1000), order.by=as.Date("2003-02-01")+1:1000)
# tfplot(z, Xaxis="auto", xlab= "time, time, time", source="Source: RNG")

 unlink("Rplots.pdf")
