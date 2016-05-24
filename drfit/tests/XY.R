library(drfit)
data(XY)
rXY <- drfit(XY,logit=TRUE,weibull=TRUE,chooseone=FALSE)
print(rXY,digits=5)
