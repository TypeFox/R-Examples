#Loading package
library(R0)

## Data is taken from the paper by Nishiura for key transmission parameters of an institutional
## outbreak during 1918 influenza pandemic in Germany)
data(Germany.1918)

## For this exemple, we use the exact same call as for the internal sensitivity analysis function

## sa.type = "GT"

## Here we will test GT with means of 1 to 5, each time with SD constant (1)
## GT and SD can be either fixed value or vectors of values
## Actual value in simulations may differ, as they are adapted according to the distribution type
tmp<-sensitivity.analysis(sa.type="GT", incid=Germany.1918, GT.type="gamma", GT.mean=seq(1,5,1), 
                          GT.sd.range=1, begin=1, end=27, est.method="EG")

## Results are stored in a matrix, each line dedicated to a (mean,sd) couple
plot(x=tmp[,"GT.Mean"], xlab="mean GT (days)", y=tmp[,"R"], ylim=c(1.2, 2.1), ylab="R0 (95% CI)", 
     type="p", pch=19, col="black", main="Sensitivity of R0 to mean GT")
arrows(x0=as.numeric(tmp[,"GT.Mean"]), y0=as.numeric(tmp[,"CI.lower"]), 
       y1=as.numeric(tmp[,"CI.upper"]), angle=90, code=3, col="black", length=0.05)
## One could tweak this example to change sorting of values (per mean, or per standard deviation)
## eg: 'x=tmp[,c('GT.Mean')]' could become 'x=tmp[,c('GT.SD')]'


## sa.type="time"

mGT<-generation.time("gamma", c(2.6,1))
sen=sensitivity.analysis(sa.type="time", incid=Germany.1918, GT=mGT, begin=1:15, end=16:30, 
                         est.method="EG")
# ...
# Warning message:
# If 'begin' and 'end' overlap, cases where begin >= end are skipped.
# These cases often return Rsquared = 1 and are thus ignored.
## A list with different estimates of reproduction ratio, exponential growth rate and 95%CI 
## wtih different pairs of begin and end dates in form of data frame is returned.
## If method is "EG", results will include growth rate and deviance R-squared measure
## Else, if "ML" method is used, growth rate and R-squared will be set as NA

## Interesting results include the variation of R0 given specific begin/end dates.
## Such results can be plot as a colored matrix and display Rsquared=f(time period)
plot(sen, what=c("criterion","heatmap"))
## Returns complete data.frame of best R0 value for each time period 
## (allows for quick visualization)
## The "best.fit" is the time period over which the estimate is the more robust

# $best.fit
#    Time.period Begin.dates  End.dates       R Growth.rate  Rsquared CI.lower. CI.upper.
# 92          15  1970-01-08 1970-01-23 1.64098   0.1478316 0.9752564  1.574953  1.710209
