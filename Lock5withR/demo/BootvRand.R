
BootvRand <- function( data, i=1, null=0, xlab="", level=.95, alpha=1-level, ... ) { 
  data$percentile <- pdata( data[,i], data[,i] )
  M <- mean( data[,i], na.rm=TRUE )
  data1 <- data
  data1$dist <- "bootstrap"
  data2 <- data
  data2$dist <- "randomization"
  data1$offset <- 0
  data2$offset <- null - M
  fancydata <- rbind(data1, data2)
  histogram( ~(fancydata[,i] + offset) | dist, 
             groups = cut(percentile, c(0,alpha/2,1-alpha/2,1)),
             data=fancydata, 
             fcol = c("red","skyblue","red"),
             v = c( mean(data1[,i], na.rm=TRUE), 
                    mean(data2[,i] + data2[,"offset"] , na.rm=TRUE)),
             layout=c(1,2), xlab=xlab,
             ...)
}

mBootvRand <- function(data, i=1, null=mean(data[,i],na.rm=TRUE), digits=2,...) {
  manipulate( BootvRand(data, i=i, null=N, alpha=ALPHA, ...),
              N = slider(
                label="Null Hypothesis Value",
                round(mean(data[,i], na.rm=TRUE) - 4 * iqr(data[,i], na.rm=TRUE), digits),
                round(mean(data[,i], na.rm=TRUE) + 4 * iqr(data[,i], na.rm=TRUE), digits),
                round(null, digits), 
                step=10^(-digits)
              ),
              ALPHA = picker(label=expression(alpha), 
							 "0.10, 90% confidence" = 0.1, 
							 "0.05, 95% confidence" = 0.05, 
							 "0.01, 99% confidence" = 0.01, 
                             initial="0.05, 95% confidence")
  )
}

require(Lock5Data)
require(mosaic)
require(manipulate)
Boot.Temp <- do(1000) * mean( ~ BodyTemp, data=resample(BodyTemp50) )
mBootvRand( Boot.Temp, null=98.6 )

