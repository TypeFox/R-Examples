## ----'preamble', include=FALSE, warning=FALSE, message=FALSE----
library(knitr)

# set the knitr options ... for everyone!
# if you unset this, then vignette build bonks. oh, joy.
#opts_knit$set(progress=TRUE)
opts_knit$set(eval.after='fig.cap')
# for a package vignette, you do want to echo.
# opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE)
opts_chunk$set(warning=FALSE,message=FALSE)
#opts_chunk$set(results="asis")
opts_chunk$set(cache=TRUE,cache.path="cache/SharpeR")

#opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
opts_chunk$set(fig.path="figure/SharpeR",dev=c("pdf"))
opts_chunk$set(fig.width=5,fig.height=4,dpi=64)

# doing this means that png files are made of figures;
# the savings is small, and it looks like shit:
#opts_chunk$set(fig.path="figure/",dev=c("png","pdf","cairo_ps"))
#opts_chunk$set(fig.width=4,fig.height=4)
# for figures? this is sweave-specific?
#opts_knit$set(eps=TRUE)

# this would be for figures:
#opts_chunk$set(out.width='.8\\textwidth')
# for text wrapping:
options(width=64,digits=2)
opts_chunk$set(size="small")
opts_chunk$set(tidy=TRUE,tidy.opts=list(width.cutoff=50,keep.blank.line=TRUE))

compile.time <- Sys.time()

# from the environment

# only recompute if FORCE_RECOMPUTE=True w/out case match.
FORCE_RECOMPUTE <- 
	(toupper(Sys.getenv('FORCE_RECOMPUTE',unset='False')) == "TRUE")

# compiler flags!

# not used yet
LONG.FORM <- FALSE

library(quantmod)
options("getSymbols.warning4.0"=FALSE)

## ----'babysteps'----------------------------------------------
library(SharpeR)
# suppose you computed the Sharpe of your strategy to
# be 1.3 / sqrt(yr), based on 1200 daily observations.
# store them as follows:
my.sr <- sr(sr=1.3,df=1200-1,ope=252,epoch="yr")
print(my.sr)
# multiple strategies can be tracked as well.
# one can attach names to them.
srstats <- c(0.5,1.2,0.6)
dim(srstats) <- c(3,1)
rownames(srstats) <- c("strat. A","strat. B","benchmark")
my.sr <- sr(srstats,df=1200-1,ope=252,epoch="yr")
print(my.sr)

## ----'showoff'------------------------------------------------
set.seed(as.integer(charToRaw("set the seed")))
# Sharpe's 'model': just given a bunch of returns.
returns <- rnorm(253*8,mean=3e-4,sd=1e-2)
asr <- as.sr(returns,ope=253,epoch="yr")
print(asr)
# a data.frame with a single strategy
asr <- as.sr(data.frame(my.strategy=returns),ope=253,epoch="yr")
print(asr)

## ----'more_data_frame'----------------------------------------
# a data.frame with multiple strategies
asr <- as.sr(data.frame(strat1=rnorm(253*8),strat2=rnorm(253*8,mean=4e-4,sd=1e-2)),
	ope=253,epoch="yr")
print(asr)

## ----'stock_loading',eval=FALSE,echo=TRUE---------------------
#  require(quantmod)
#  # get price data, compute log returns on adjusted closes
#  get.ret <- function(sym,warnings=FALSE,...) {
#  	# getSymbols.yahoo will barf sometimes; do a trycatch
#    trynum <- 0
#  	while (!exists("OHCLV") && (trynum < 7)) {
#  		trynum <- trynum + 1
#  		try(OHLCV <- getSymbols(sym,auto.assign=FALSE,warnings=warnings,...),silent=TRUE)
#    }
#  	adj.names <- paste(c(sym,"Adjusted"),collapse=".",sep="")
#  	if (adj.names %in% colnames(OHLCV)) {
#  		adj.close <- OHLCV[,adj.names]
#  	} else {
#  		# for DJIA from FRED, say.
#  		adj.close <- OHLCV[,sym]
#  	}
#  	rm(OHLCV)
#  	# rename it
#  	colnames(adj.close) <- c(sym)
#  	adj.close <- adj.close[!is.na(adj.close)]
#  	lrets <- diff(log(adj.close))
#  	#chop first
#  	lrets[-1,]
#  }
#  get.rets <- function(syms,...) { some.rets <- do.call("cbind",lapply(syms,get.ret,...)) }

## ----'stock_loading_sneaky',eval=TRUE,echo=FALSE--------------
# sleight of hand to load precomputed data instead.
get.rets <- function(syms,from='2003-01-01',to='2013-01-01',...) {
	fname <- system.file('extdata','ret_data.rda',package='SharpeR')
	if (fname == "") {
		fname <- 'ret_data.rda'
	}
	# poofs all.ret here
	load(fname)
	sub.data <- all.ret[paste(from,to,sep="::"),colnames(all.ret) %in% syms]
	return(sub.data)
}

## ----'helper_function'----------------------------------------
require(quantmod)
# quantmod::periodReturn does not deal properly with multiple
# columns, and the straightforward apply(mtms,2,periodReturn) barfs
my.periodReturn <- function(mtms,...) {
	per.rets <- do.call(cbind,lapply(mtms,
		function(x) {
			retv <- periodReturn(x,...)
			colnames(retv) <- colnames(x)
			return(retv) 
		}))
}
# convert log return to mtm, ignoring NA
lr2mtm <- function(x,...) {
	x[is.na(x)] = 0
	exp(cumsum(x))
}

## ----'some_stocks'--------------------------------------------
some.rets <- get.rets(c("IBM","AAPL","XOM"),
	from="2007-01-01",to="2013-01-01")
print(as.sr(some.rets))

## ----'reannualize'--------------------------------------------
yearly <- as.sr(some.rets[,"XOM"])
monthly <- reannualize(yearly,new.ope=21,new.epoch="mo.")
print(yearly)
# significance should be the same, but units changed.
print(monthly)

## ----'AAPL'---------------------------------------------------
# get the returns (see above for the function)
aapl.rets <- get.rets(c("AAPL","SPY"),from="2003-01-01",to="2013-01-01")
# make them monthly:
mo.rets <- my.periodReturn(lr2mtm(aapl.rets),period='monthly',type='arithmetic')
rm(aapl.rets)  # cleanup
# look at both of them together:
both.sr <- as.sr(mo.rets)
print(both.sr)
# confindence intervals on the Sharpe:
print(confint(both.sr))
# perform a CAPM attribution, using SPY as 'the market'
linmod <- lm(AAPL ~ SPY,data=mo.rets)
# convert attribution model to Sharpe
CAPM.sr <- as.sr(linmod,ope=both.sr$ope,epoch="yr")
# statistical significance does not change (though note the sr summary
# prints a 1-sided p-value)
print(summary(linmod))
print(CAPM.sr)
# the confidence intervals tell the same story, but in different units:
print(confint(linmod,'(Intercept)'))
print(confint(CAPM.sr))

## ----'SPDRcheck'----------------------------------------------
# get the sector 'spiders'
secto.rets <- get.rets(c("XLY","XLE","XLP","XLF","XLV","XLI","XLB","XLK","XLU"),
	from="2003-01-01",to="2013-01-01")
# make them monthly:
mo.rets <- my.periodReturn(lr2mtm(secto.rets),period='monthly',type='arithmetic')
# one-sample test on utilities:
XLU.monthly <- mo.rets[,"XLU"]
print(sr_test(XLU.monthly),alternative="two.sided")

# test for equality of Sharpe among the different spiders
print(sr_equality_test(secto.rets))
# perform a paired two-sample test via sr_test:
XLF.monthly <- mo.rets[,"XLF"]
print(sr_test(x=XLU.monthly,y=XLF.monthly,ope=12,paired=TRUE))

## ----'sropt_basics'-------------------------------------------
set.seed(as.integer(charToRaw("7bf4b86a-1834-4b58-9eff-6c7dec724fec")))
# from a matrix object:
ope <- 253
n.stok <- 7
n.yr <- 8
# somewhat unrealistic: independent returns.
rand.rets <- matrix(rnorm(n.yr * ope * n.stok),ncol=n.stok)
asro <- as.sropt(rand.rets,ope=ope)
rm(rand.rets)
print(asro)
# under the alternative, when the mean is nonzero
rand.rets <- matrix(rnorm(n.yr * ope * n.stok,mean=6e-4,sd=1e-2),ncol=n.stok)
asro <- as.sropt(rand.rets,ope=ope)
rm(rand.rets)
print(asro)
# from an xts object
some.rets <- get.rets(c("IBM","AAPL","XOM"),
	from="2007-01-01",to="2013-01-01")
asro <- as.sropt(some.rets)
print(asro)

## ----'sropt_estim'--------------------------------------------
# confidence intervals:
print(confint(asro,level.lo=0.05,level.hi=1))
# estimation
print(inference(asro,type="KRS"))
print(inference(asro,type="MLE"))

## ----'MLE_rule'-----------------------------------------------
ope <- 253
zeta.s <- 0.8
n.check <- 1000
df1 <- 10
df2 <- 6 * ope
rvs <- rsropt(n.check,df1,df2,zeta.s,ope,drag=0)
roll.own <- sropt(z.s=rvs,df1,df2,drag=0,ope=ope,epoch="yr")
MLEs <- inference(roll.own,type="MLE")
zerMLE <- MLEs <= 0
crit.value <- 0.5 * (max(rvs[zerMLE])^2 + min(rvs[!zerMLE])^2)
aspect.ratio <- df1 / (df2 / ope)
cat(sprintf("empirical cutoff for zero MLE is %2.2f yr^{-1}\n", crit.value))
cat(sprintf("the aspect ratio is %2.2f yr^{-1}\n",aspect.ratio))

## ----'estimate_overfit',fig.cap=paste("Q-Q plot of",n.sim,"achieved optimal \\txtSR values from brute force search over both windows of a Moving Average Crossover under the null of driftless log returns with zero autocorrelation versus the approximation by a 2-parameter optimal \\txtSR distribution is shown.")----
require(TTR)
# brute force search two window MAC
brute.force <- function(lrets,rrets=exp(lrets)-1,win1,win2=win1) {
	mtms <- c(1,exp(cumsum(lrets)))  # prepend a 1.
  # do all the SMAs;
  SMA1 <- sapply(win1,function(n) { SMA(mtms,n=n) }) 
  symmetric <- missing(win2)
  if (!symmetric)
  	SMA2 <- sapply(win2,function(n) { SMA(mtms,n=n) }) 

  mwin <- max(c(win1,win2))
  zeds <- matrix(NaN,nrow=length(win1),ncol=length(win2))
	upb <- if (symmetric) length(win1) - 1 else length(win1)
	# 2FIX: vectorize this!
	for (iidx in 1:upb) {
		SM1 <- SMA1[,iidx]
		lob <- if (symmetric) iidx + 1 else 1
		for (jidx in lob:length(win2)) {
			SM2 <- if (symmetric) SMA1[,jidx] else SMA2[,jidx]
			trades <- sign(SM1 - SM2)
			dum.bt <- trades[mwin:(length(trades)-1)] * rrets[mwin:length(rrets)]  # braindead backtest.
			mysr <- as.sr(dum.bt)
			zeds[iidx,jidx] <- mysr$sr
			# abuse symmetry of arithmetic returns
			if (symmetric) zeds[jidx,iidx] <- - zeds[iidx,jidx]  
		}
	}
	retv <- max(zeds,na.rm=TRUE) 
	return(retv)
}
# simulate one.
sim.one <- function(nbt,win1,...) {
	lrets <- rnorm(nbt+max(win1),sd=0.01)
	retv <- brute.force(lrets,win1=win1,...)
	return(retv)
}
# set everything up
set.seed(as.integer(charToRaw("e23769f4-94f8-4c36-bca1-28c48c49b4fb")))
ope <- 253
n.yr <- 4
n.obs <- ceiling(ope * n.yr)
LONG.FORM <- FALSE
n.sim <- if (LONG.FORM) 2048 else 1024
win1 <- if (LONG.FORM) c(2,4,8,16,32,64,128,256) else c(4,16,64,256)

# run them
system.time(max.zeds <- replicate(n.sim,sim.one(n.obs,win1)))
# qqplot;
qqplot(qsropt(ppoints(length(max.zeds)),df1=2,df2=n.obs),max.zeds,
			 xlab = "Theoretical Approximate Quantiles", ylab = "Sample Quantiles")
qqline(max.zeds,datax=FALSE,distribution = function(p) { qsropt(p,df1=2,df2=n.obs) },
			 col=2)

## ----'now_on_spy'---------------------------------------------
# is MAC on SPY significant?
SPY.lret <- get.rets(c('SPY'),from="2003-01-01",to="2013-01-01")
# oops! there might be NAs in there!
mysr <- as.sr(SPY.lret,na.rm=TRUE)  # just to get the ope
print(mysr)
# get rid of NAs
SPY.lret[is.na(SPY.lret)] <- 0
# try a whole lot of windows:
win1 <- seq(4,204,by=10)
zeds <- brute.force(SPY.lret,win1=win1)
SPY.MAC.asro <- sropt(z.s=zeds,df1=2,df2=length(SPY.lret) - max(win1),ope=mysr$ope)
print(SPY.MAC.asro)
print(inference(SPY.MAC.asro,type="KRS"))

## ----'del_sropt_basics'---------------------------------------
set.seed(as.integer(charToRaw("364e72ab-1570-43bf-a1c6-ee7481e1c631")))
# from a matrix object:
ope <- 253
n.stok <- 7
n.yr <- 8
# somewhat unrealistic: independent returns, under the null
rand.rets <- matrix(rnorm(n.yr * ope * n.stok),ncol=n.stok)
# the hedge constraint: hedge out the first stock.
G <- diag(n.stok)[1,]
asro <- as.del_sropt(rand.rets,G,ope=ope)
print(asro)
# hedge out the first two 
G <- diag(n.stok)[1:2,]
asro <- as.del_sropt(rand.rets,G,ope=ope)
print(asro)
# under the alternative, when the mean is nonzero
rand.rets <- matrix(rnorm(n.yr * ope * n.stok,mean=6e-4,sd=1e-2),ncol=n.stok)
G <- diag(n.stok)[1,]
asro <- as.del_sropt(rand.rets,G,ope=ope)
print(asro)

## ----'del_sropt_hedging'--------------------------------------
# from an xts object
some.rets <- get.rets(c("SPY","IBM","AAPL","XOM"),
	from="2007-01-01",to="2013-01-01")
# without the hedge, allowing SPY position
asro <- as.sropt(some.rets)
print(asro)
# hedge out SPY!
G <- diag(dim(some.rets)[2])[1,]
asro.hej <- as.del_sropt(some.rets,G)
print(asro.hej)

## ----'del_sropt_estim'----------------------------------------
# confidence intervals:
print(confint(asro,level.lo=0.05,level.hi=1))
# estimation
print(inference(asro,type="KRS"))
print(inference(asro,type="MLE"))

## ----'sr_eq_test_normal',fig.cap=paste("P-values from \\Rfunction{sr\\_equality\\_test} on",n.sim,"multiple realizations of independent, normally distributed returns series are Q-Q plotted against uniformity.")----
# function for Q-Q plot against uniformity
qqunif <- function(x,xlab="Theoretical Quantiles under Uniformity",
									 ylab=NULL,...) {
	if (is.null(ylab))
		ylab=paste("Sample Quantiles (",deparse(substitute(x)),")",sep="")
	qqplot(qunif(ppoints(length(x))),x,xlab=xlab,ylab=ylab,...)
	abline(0,1,col='red')
}

# under normality.
LONG.FORM <- FALSE
n.ro <- if (LONG.FORM) 253*4 else 253*2
n.co <- if (LONG.FORM) 20 else 4
n.sim <- if (LONG.FORM) 1024 else 512
set.seed(as.integer(charToRaw("e3709e11-37e0-449b-bcdf-9271fb1666e5")))
afoo <- replicate(n.sim,sr_equality_test(matrix(rnorm(n.ro*n.co),ncol=n.co),type="F"))
sr_eq.pvals <- unlist(afoo["p.value",])
qqunif(sr_eq.pvals)

## ----'sr_eq_test_hetero',fig.cap=paste("P-values from \\Rfunction{sr\\_equality\\_test} on",n.sim,"multiple realizations of independent, \\tlaw{}-distributed returns series are Q-Q plotted against uniformity.")----
# try heteroskedasticity?
set.seed(as.integer(charToRaw("81c97c5e-7b21-4672-8140-bd01d98d1d2e")))
afoo <- replicate(n.sim,sr_equality_test(matrix(rt(n.ro*n.co,df=4),ncol=n.co),type="F"))
sr_eq.pvals <- unlist(afoo["p.value",])
qqunif(sr_eq.pvals)

## ----'sr_eq_test_correlated',fig.cap=paste("P-values from \\Rfunction{sr\\_equality\\_test} on",n.sim,"multiple realizations of correlated, normally-distributed returns series are Q-Q plotted against uniformity.")----
# try correlated returns
n.fact <- max(2,n.co - 5)
gen.ret <- function(n1,n2,f=max(2,n2-2),fuzz=0.1) {
	A <- matrix(rnorm(n1*f),nrow=n1)
	B <- matrix(rnorm(f*n2),nrow=f)
	C <- sqrt(1-fuzz^2) * A %*% B + fuzz * matrix(rnorm(n1*n2),nrow=n1)
}
set.seed(as.integer(charToRaw("e4d61c2c-efb3-4cba-9a6e-5f5276ce2ded")))
afoo <- replicate(n.sim,sr_equality_test(gen.ret(n.ro,n.co,n.fact),type="F"))
sr_eq.pvals <- unlist(afoo["p.value",])
qqunif(sr_eq.pvals)

## ----'sr_eq_test_mtime_pre',echo=FALSE------------------------
n.co <- if (LONG.FORM) 100 else 50 
n.sim <- if (LONG.FORM) 1024 else 128

## ----'sr_eq_test_mtime',fig.cap=paste("P-values from \\Rfunction{sr\\_equality\\_test} on",n.sim,"multiple realizations of",n.co,"market timing strategies' returns series are Q-Q plotted against uniformity. Nominal coverage is not maintained: the test is far too liberal in this case.")----
n.co <- if (LONG.FORM) 100 else 50 
n.sim <- if (LONG.FORM) 1024 else 128
SPY.lret <- get.rets(c('SPY'),from="2003-01-01",to="2013-01-01")
SPY.wk.rret <- my.periodReturn(lr2mtm(SPY.lret),period='weekly',
	type='arithmetic')
gen.tim <- function(n2) {
	mkt.timing.signal <- sign(rnorm(n2*length(SPY.wk.rret)))
	mkt.ret <- matrix(rep(SPY.wk.rret,n2) * mkt.timing.signal,ncol=n2)
}
set.seed(as.integer(charToRaw("447cfe85-b612-4b14-bd01-404e6e99aca4")))
system.time(afoo <- replicate(n.sim,sr_equality_test(gen.tim(n.co),type="F")))
sr_eq.pvals <- unlist(afoo["p.value",])
qqunif(sr_eq.pvals)

## ----'sr_fancy_eq'--------------------------------------------
# get returns
some.rets <- get.rets(c("IBM","AAPL","XOM"),
	from="2007-01-01",to="2013-01-01")
# using the default vcov
test.vanilla <- sr_equality_test(some.rets,type="F")
print(test.vanilla)
if (require(sandwich)) {
	# and a fancy one:
	test.HAC <- sr_equality_test(some.rets,type="F",vcov.func=vcovHAC)
	print(test.HAC)
}

## ----'sr_vcov'------------------------------------------------
# get returns
some.rets <- get.rets(c("IBM","AAPL","XOM"),
	from="2007-01-01",to="2013-01-01")
ope <- 252
vanilla.Sig <- sr_vcov(some.rets,ope=ope)
print(vanilla.Sig)
if (require(sandwich)) {
	HC.Sig <- sr_vcov(some.rets,vcov=vcovHC,ope=ope)
	print(HC.Sig$Ohat)
}
if (require(sandwich)) {
	HAC.Sig <- sr_vcov(some.rets,vcov=vcovHAC,ope=ope)
	print(HAC.Sig$Ohat)
}

## ----'marko_vcov'---------------------------------------------
# get returns
some.rets <- get.rets(c("IBM","AAPL","XOM"),
	from="2007-01-01",to="2013-01-01")

ism.wald <- function(X,vcov.func=vcov) {
	# negating returns is idiomatic to get + Markowitz
	ism <- ism_vcov(- as.matrix(X),vcov.func=vcov.func)
	ism.mu <- ism$mu[1:ism$p]
	ism.Sg <- ism$Ohat[1:ism$p,1:ism$p]
	retval <- ism.mu / sqrt(diag(ism.Sg))
	dim(retval) <- c(ism$p,1)
	rownames(retval) <- rownames(ism$mu)[1:ism$p]
	return(retval)
}

wald.stats <- ism.wald(some.rets)
print(t(wald.stats))

if (require(sandwich)) {
	wald.stats <- ism.wald(some.rets,vcov.func=sandwich::vcovHAC)
	print(t(wald.stats))
}

