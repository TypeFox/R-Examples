#rm(list=ls())
#source("~/lavori/Rdevel/SGR1.0/R/dgBetaD.R")
#k <- 7
#h <- 7
#Q <- 7
#gam <- 1; del <- 1; pi <- .5

pfakebad <-
function(k,h=k,p=0,Q=5,gam=1,del=1) {
  if ((k>Q)|(h>Q)) {
    warning("value/s out of scale boundaries")
    outP <- 0
  }
  gB <- dgBetaD(k,1,h-1,gam,del)
  if ((h==k)&(h==1)) {
    outP <- 1
  } else {
    if ((1<=k)&(k<h)&(h<=Q)) outP <- p*gB
    if ((1<k)&(k==h)&(h<=Q)) outP <- 1-p
    if ((1<=h)&(h<k)&(k<=Q)) outP <- 0
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
#  for (i in x) y <- c(y,pfakebad(x[i],h=5,Q=7,gam=GA[j],del=DE[j],p=.4))
#  plot(x,y,type="h",panel.first=points(x,y,pch=19),
#       main=paste("gamma=",GA[j]," delta=",DE[j],sep=""),ylim=c(0,.7),
#       ylab="Replacement probability")  
#}