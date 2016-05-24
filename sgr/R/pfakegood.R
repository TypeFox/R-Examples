#rm(list=ls())
#source("~/lavori/Rdevel/SGR1.0/R/dgBetaD.R")

pfakegood <-
function(k,h=k,p=0,Q=5,gam=1,del=1) {
  if ((k>Q)|(h>Q)) {
    warning("value/s out of scale boundaries")
    outP <- 0
  }
  gB <- dgBetaD(k,h+1,Q,gam,del)
  if ((h==k)&(h==Q)) {
    outP <- 1
  } else {
    if ((1<=h)&(h<k)&(k<=Q)) outP <- p*gB
    if ((1<=k)&(k==h)&(h<Q)) outP <- 1-p
    if ((1<=k)&(k<h)&(h<=Q)) outP <- 0
  }
  return(outP)
  #return(list(outP=outP,gBetaD=gB,Gd=gB1))
}

## example
#x <- 1:7
#GA <- c(1,3,1.5,8); DE <- c(1,3,4,2.5)
#par(mfrow=c(2,2))
#for (j in 1:4) {
#  y <- NULL
#  for (i in x) y <- c(y,pfakegood(x[i],h=3,Q=7,gam=GA[j],del=DE[j],alpha=.4))
#  plot(x,y,type="h",panel.first=points(x,y,pch=19),
#       main=paste("gamma=",GA[j]," delta=",DE[j],sep=""),ylim=c(0,.7),
#       ylab="Replacement probability")  
#}