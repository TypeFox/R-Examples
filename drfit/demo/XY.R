require(drfit)
data(XY)
rXY <- drfit(XY)
rXY
drplot(rXY,XY,dtype="raw",
    overlay=TRUE,bw=FALSE)
