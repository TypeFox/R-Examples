library(RobPer)

################
## betaCvMfit ##
################

rm(list=ls())

n <- 100
set.seed(1066)
dat1 <- rbeta(n, 3, 2)
dat2 <- dat1
set.seed(1216)
ii <- sample(1:n, n/10)
set.seed(1399)
dat2[ii] <- rbeta(n/10, 1, 4)

betaCvMfit(dat1, CvM=TRUE, rob=TRUE)
betaCvMfit(dat1, CvM=TRUE, rob=FALSE)
betaCvMfit(dat1, CvM=FALSE, rob=TRUE)
betaCvMfit(dat1, CvM=FALSE, rob=FALSE)

betaCvMfit(dat2, CvM=TRUE, rob=TRUE)
betaCvMfit(dat2, CvM=TRUE, rob=FALSE)
betaCvMfit(dat2, CvM=FALSE, rob=TRUE)
betaCvMfit(dat2, CvM=FALSE, rob=FALSE)


############
## RobPer ##
############

## (with application of auxiliary functions checkbest, FastS, FastTau, IRWLS, IWLSiteration, re.s, singlePerM, singlePernotM and Xgen)

rm(list=ls())

set.seed(1066)
## (also application of tsgen and some of its auxiliary functions)
ts <- tsgen(ttype="sine", ytype="trian", pf=7, redpart=0.1, s.outlier.fraction=0.1, interval=TRUE, npoints=100, ncycles=5, ps=10, SNR=0.5)

per <- 6:15

set.seed(1216)
ww <- sample(c(TRUE, FALSE), 13, TRUE)

set.seed(1485)
seeds <- sample(1:1399, 13, TRUE)
i <- 1

for (reg in c("L2", "L1", "huber", "bisquare", "S", "tau")) {
    set.seed(seeds[i])
    print(signif(RobPer(ts, weighting=ww[i], periods=per, regression=reg, model="sine"), 3))
    i <- i+1
}

set.seed(seeds[i])
signif(RobPer(ts, weighting=ww[i], periods=per, regression="LTS", model="sine", LTSopt=FALSE), 3)
i <- i+1
set.seed(seeds[i])
signif(RobPer(ts, weighting=ww[i], periods=per, regression="LTS", model="sine", LTSopt=TRUE), 3)
i <- i+1

for (mod in c("step", "2step", "fourier(2)", "fourier(3)", "splines")) {
    set.seed(seeds[i])
    print(signif(RobPer(ts, weighting=ww[i], periods=per, regression="huber", model=mod, steps=5), 3))
    i <- i+1
}


###########
## tsgen ##
###########

## (with application of auxiliary functions disturber, lc_noise, sampler, signalgen, TK95 and TK95_uneq)

rm(list=ls())

pf <- 7
rp <- 0.1
sof <- 0.1
n <- 100
nc <- 5
ps <- 20
snr <- 0.5

set.seed(1066)
seeds <- matrix(sample(1:1216, 32, T), nrow=16)
i <- 1

for (tt in c("equi", "unif", "sine", "trian")) for (yt in c("const", "sine", "trian", "peak")) {
    set.seed(seeds[i,1])
    print(tsgen(ttype=tt, ytype=yt, pf=pf, redpart=rp, s.outlier.fraction=sof, interval=TRUE, npoints=n, ncycles=nc, ps=ps, SNR=snr))
    set.seed(seeds[i,2])
    print(tsgen(ttype=tt, ytype=yt, pf=pf, redpart=rp, s.outlier.fraction=sof, interval=FALSE, npoints=n, ncycles=nc, ps=ps, SNR=snr))
    i <- i+1
}