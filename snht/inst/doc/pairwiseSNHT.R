## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----setup, include=FALSE, cache=FALSE, echo=FALSE----------------------------
library(knitr)
opts_chunk$set(fig.path = 'figure/', fig.align = 'center', fig.show = 'hold',
               warning = FALSE, message = FALSE, error = FALSE, tidy = FALSE,
               results = 'markup', eval = TRUE, echo = TRUE, cache = FALSE)
options(replace.assign = TRUE, width = 80)

## ----eval=FALSE---------------------------------------------------------------
#  pairwiseSNHT(data, dist, k, period, crit, returnStat=TRUE/FALSE)

## -----------------------------------------------------------------------------
library(snht)
library(reshape2)
library(ggplot2)
library(mvtnorm)

## -----------------------------------------------------------------------------
set.seed(2)
Cor<-rbind(c(0.5,0.8,0.8,0.8,0.8),c(0,0.5,0.8,0.5,0.6),
c(0,0,0.5,0.8,0.5),c(0,0,0,0.5,0.6),c(0,0,0,0,0.5))
Cor<-t(Cor)+Cor
baseData<-rmvnorm(mean=rep(0,5),sig=Cor,n=1000)+cos(1:1000*2*pi/200)
baseData[401:1000,1]<-baseData[401:1000,1]+0.5

## ----echo = FALSE-------------------------------------------------------------
toPlot = data.frame(baseData)
toPlot$Time = 1:1000
toPlot = melt(toPlot, id.vars = "Time")
toPlot$variable = gsub("X", "Series ", toPlot$variable)
ggplot(toPlot, aes(x = Time)) +
    geom_line(aes(y = value, color = variable)) +
    facet_grid(variable ~ .) + guides(color = FALSE)

## -----------------------------------------------------------------------------
dist<-matrix(0,5,5)
dist<-dist(rbind(c(1,0),c(1,1),c(0,0),c(1,-1),c(2,0)))
dist<-as.matrix(dist)
colnames(dist)<-rownames(dist)<-1:5
dist

## -----------------------------------------------------------------------------
colnames(baseData) <-"1":"5"
baseData <- data.frame(time = 1:1000, baseData)
baseData <- melt(baseData, id.vars = "time", variable.name = "location",
                 value.name = "data")
baseData$location<-gsub("X","",baseData$location)

out1 <- pairwiseSNHT(baseData, dist, k=3, period=200, 
crit=qchisq(1-0.05/600,df=1), returnStat=T)

pairs <- colnames(out1)
pairs
out1[300:302, ]

## ----echo = FALSE-------------------------------------------------------------
par(mfrow=c(4,2),mar=c(1,1,1,1))
for(i in 1:8){
  plot(out1[,i],type="l",ylim=c(0,100),
       xlab="time",ylab="SNHT statistic")
  title(pairs[i], line = -1) 
  abline(h=qchisq(1-0.05/600,df=1), col="red")
}

## -----------------------------------------------------------------------------
out2 <- pairwiseSNHT(baseData, dist, k=3, period=200, 
crit=qchisq(1-0.05/600,df=1), returnStat=F)
out2$breaks
str(out2$data)

## -----------------------------------------------------------------------------
newPair2 <- out2$data
outNew1 <- pairwiseSNHT(newPair2,dist,k=3,period=200,
                        crit=qchisq(1-0.05/600,df=1),returnStat=T)

## ----echo = FALSE-------------------------------------------------------------
par(mfrow=c(4,2),mar=c(1,1,1,1))
for(i in 1:8){
  plot(outNew1[,i],type="l",ylim=c(0,100),
       xlab="time",ylab="SNHT statistic")
  title(pairs[i], line = -1) 
  abline(h=qchisq(1-0.05/600,df=1), col="red")
}

