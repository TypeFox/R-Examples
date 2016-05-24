require(drfit)
data(IM1xVibrio)
rIM1xVibrio <- drfit(IM1xVibrio)
rIM1xVibrio
drplot(rIM1xVibrio,IM1xVibrio,
    overlay=TRUE,bw=FALSE,lpos="bottomleft")
