#rm(list=ls())
#source("~/lavori/Rdevel/sgr1.0/R/dgBetaD.R")
#source('~/lavori/Rdevel/sgr1.0/R/model.fake.par.R')
#k <- 6
#h <- 7
#Q <- 7
#gam <- c(1,1); del <- c(1,1); p <- c(.4,0)
#fakemodel=c(NA,"uninformative","average","slight","extreme")

pfake <-
function(k,h=k,p=c(0,0),Q=5,gam=c(1,1),del=c(1,1),
         fake.model=c("uninformative","average","slight","extreme")) {
  
  fake.model <- match.arg(fake.model)
  if (fake.model!="uninformative") {
    MF <- model.fake.par(fake.model)
    gam <- MF$gam; del <- MF$del
  }
  
  if (sum(p)>1) {
    warning("sum(p) cannot be greater than 1")
    p <- p/sum(p) 
  }
  
  if ((k>Q)|(h>Q)) {
    warning("value/s out of scale boundaries")
    outP <- 0
  }

  if ((1<=k)&(k<h)&(h<=Q)) outP <- p[2]*dgBetaD(k,1,h-1,gam[2],del[2])  
  if ((1<=k)&(k==h)&(h<=Q)) {
    outP <- 1-(p[1]+p[2])
    if (k==1) outP <- 1-p[1]
    if (k==Q) outP <- 1-p[2]
  }
  if ((1<=h)&(h<k)&(k<=Q)) outP <- p[1]*dgBetaD(k,h+1,Q,gam[1],del[1])

  return(outP)
}

## example
#x <- 1:7
#GA <- c(1,3,1.5,8); DE <- c(1,3,4,2.5)
#P <- c(0.4,.4); H <- 7
#par(mfrow=c(2,2))
#for (j in 1:4) {
#  y <- NULL
#  for (i in x) y <- c(y,pfake(x[i],h=H,Q=max(x),
#                gam=c(GA[j],GA[j]),del=c(DE[j],DE[j]),p=P))
#  plot(x,y,type="h",panel.first=points(x,y,pch=19),
#       main=paste("gamma=",GA[j]," delta=",DE[j],sep=""),ylim=c(0,1),
#       ylab="Replacement probability") 
#  print(sum(y))
#}

#x <- 1:5 
#models <- c("uninformative","average","slight","extreme")
#par(mfrow=c(2,2))
#for (j in 1:4) {
#  y <- NULL
#  for (i in x) y <- c(y,pfake(x[i],h=2,Q=max(x),
#            fake.model=models[j],p=c(.45,0)))
#  plot(x,y,type="h",panel.first=points(x,y,pch=19),
#       main=paste(models[j],"model"),ylim=c(0,1),
#       ylab="Replacement probability") 
#}

#par(mfrow=c(2,2))
#for (j in 1:4) {
#  y <- NULL
#  for (i in x) y <- c(y,pfake(x[i],h=4,Q=max(x),
#                              fake.model=models[j],p=c(0,.45)))
#  plot(x,y,type="h",panel.first=points(x,y,pch=19),
#       main=paste(models[j],"model"),ylim=c(0,1),
#       ylab="Replacement probability") 
#}

## nota: se sei all'estremo, h=7 che senso ha che uno dei due fake sia diverso da 0??
## quindi se sono valori estremi deve azzerarne uno
#P = c(0,.4)
#par(mfrow=c(2,4))
#for (j in x) {
#  y <- NULL
#  for (i in x) {
#    y <- c(y,pfake(x[i],h=x[j],Q=max(x),gam=c(GA[1],GA[1]),del=c(DE[1],DE[1]),p=P))
#  }
#  plot(x,y,type="h",panel.first=points(x,y,pch=19),
#       main=paste("h=",x[j],sep=""),ylim=c(0,1),
#       ylab="Replacement probability") 
#  print(sum(y,na.rm=TRUE)) 
#}