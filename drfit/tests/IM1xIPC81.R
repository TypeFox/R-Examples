library(drfit)
data(IM1xIPC81)
rIM1xIPC81 <- drfit(IM1xIPC81,linlogit=TRUE)
print(rIM1xIPC81,digits=4)
