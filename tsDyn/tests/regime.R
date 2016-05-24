library(tsDyn)

mod.set <- setar(lynx,m=2)

regime(mod.set)
regime(mod.set, initVal=FALSE)
regime(mod.set, time=FALSE)
regime(mod.set, time=FALSE, initVal=FALSE)
plot(regime(mod.set))

data(barry)
mod.tv <- TVAR(barry[,1:2], lag=2, nthresh=2, thDelay=1, trim=0.1, mTh=1, plot=TRUE)

regime(mod.tv)
regime(mod.tv, initVal=FALSE)
regime(mod.tv, time=FALSE)
regime(mod.tv, time=FALSE, initVal=FALSE)
plot(regime(mod.tv))



mod.tvecm <-TVECM(barry[,1:2], nthresh=2,lag=1, ngridBeta=20, ngridTh=30,trim=0.05, common="All")

regime(mod.tvecm)
regime(mod.tvecm, initVal=FALSE)
regime(mod.tvecm, time=FALSE)
regime(mod.tvecm, time=FALSE, initVal=FALSE)
plot(regime(mod.tvecm))
