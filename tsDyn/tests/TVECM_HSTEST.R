library(tsDyn)
data(zeroyld)
data<-zeroyld

## Test against paper:
all.equal(round(TVECM.HStest(data, lag=1, intercept=TRUE, nboot=0)$stat,4),20.5994)
all.equal(round(TVECM.HStest(data, lag=2, intercept=TRUE, nboot=0)$stat,4),28.2562 )
all.equal(round(TVECM.HStest(data, lag=3, intercept=TRUE, nboot=0)$stat,4), 29.9405 )


## prob:
all.equal(round(TVECM.HStest(data, lag=2, intercept=TRUE, nboot=0, fixed.beta=1)$stat,4),29.5295)
all.equal(round(TVECM.HStest(data, lag=1, intercept=TRUE, nboot=0, fixed.beta=1)$stat,4),21.5586 )
  
## Test: no boot
TVECM.HStest(data, lag=1, intercept=TRUE, ngridTh=50, nboot=0)
TVECM.HStest(data, lag=1, intercept=FALSE, ngridTh=50, nboot=0)
TVECM.HStest(data, lag=1, intercept=TRUE, nboot=0)
TVECM.HStest(data, lag=1, intercept=FALSE, nboot=0)


## Test: boot
set.seed(123)
t1<-TVECM.HStest(data, lag=1, intercept=TRUE, ngridTh=50, nboot=5)
set.seed(123)
t2<-TVECM.HStest(data, lag=1, intercept=FALSE, ngridTh=50, nboot=5)
set.seed(123)
t3<-TVECM.HStest(data, lag=1, intercept=TRUE, ngridTh=50, nboot=5, boot.type="ResBoot")
set.seed(123)
t4<-TVECM.HStest(data, lag=1, intercept=FALSE, ngridTh=50, nboot=5, boot.type="ResBoot")

## Test: methodst1
summary(t1)
plot(t1)
plot(t1, which="Density")
plot(t1, which="LM values")
t2
summary(t2)
t3
summary(t3)
t4
summary(t4)
