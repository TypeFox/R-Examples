source("normGibbs.r")
y<-scan("Example15.csv")
#debug(normGibbs)
normGibbs(y, priorMu = c(20, 1), priorVar = c(11.37, 1))
