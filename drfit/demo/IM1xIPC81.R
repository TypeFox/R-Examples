require(drfit)
data(IM1xIPC81)
rIM1xIPC81 <- drfit(IM1xIPC81,linlogit=TRUE)
rIM1xIPC81
drplot(rIM1xIPC81,IM1xIPC81,
    overlay=TRUE,bw=FALSE,xlim=c("auto",5.5))
