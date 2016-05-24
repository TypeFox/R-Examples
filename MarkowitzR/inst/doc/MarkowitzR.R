## ----'preamble', include=FALSE, warning=FALSE, message=FALSE, cache=FALSE----
library(knitr)

# set the knitr options ... for everyone!
# if you unset this, then vignette build bonks. oh, joy.
#opts_knit$set(progress=TRUE)
opts_knit$set(eval.after='fig.cap')
# for a package vignette, you do want to echo.
# opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE)
opts_chunk$set(warning=FALSE,message=FALSE)
#opts_chunk$set(results="asis")
opts_chunk$set(cache=TRUE,cache.path="cache/MarkowitzR")

#opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
opts_chunk$set(fig.path="figure/MarkowitzR",dev=c("pdf"))
opts_chunk$set(fig.width=4.5,fig.height=4,dpi=300)

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

mc.resolution <- ifelse(LONG.FORM,1000,200)
mc.resolution <- max(mc.resolution,100)

library(MarkowitzR)

gen_norm <- rnorm
lseq <- function(from,to,length.out) { 
	exp(seq(log(from),log(to),length.out = length.out))
}

## ----'the_proc',echo=TRUE,cache=FALSE-------------------------
require(gtools)
# set the colnames of X appropriately
set.coln <- defmacro(X,expr={
	if (is.null(colnames(X))) {
		colnames(X) <- paste0(deparse(substitute(X),nlines=1),1:(dim(X)[2]))
  }})
# compute Wald scores via this trick
wrap.itheta <- function(rets,ret.what=c('wald','mp'),...) {
	set.coln(rets)
	ret.what <- match.arg(ret.what)
	# 'flipping the sign on returns is idiomatic'
	asymv <- itheta_vcov(- as.matrix(rets),...)
	qidx <- 2:asymv$pp
	mp <- asymv$mu[qidx] 
  stat <- switch(ret.what,
		mp = mp,
		wald = mp / sqrt(diag(asymv$Ohat[2:asymv$pp,2:asymv$pp])))
  names(stat) <- colnames(rets)
	return(stat)
}
wrap.mp <- function(rets,ret.what=c('wald','mp'),...) {
	set.coln(rets)
	ret.what <- match.arg(ret.what)
	asymv <- mp_vcov(as.matrix(rets),...)
	mp <- t(asymv$W)
  stat <- switch(ret.what,
		mp = mp,
		wald = mp / sqrt(diag(asymv$What)))
	return(stat)
}
# t-stat via Britten-Jones procedure
bjones.ts <- function(rets) {
	set.coln(rets)
	ones.vec <- matrix(1,nrow=dim(rets)[1],ncol=1)
	rets <- as.matrix(rets)
	bjones.mod <- lm(ones.vec ~ rets - 1)
	bjones.sum <- summary(bjones.mod)
	retval <- bjones.sum$coefficients[,3]
  names(retval) <- colnames(rets)
	return(retval)
}
# compare the procedures
do.both <- function(rets,...) {
	set.coln(rets)
	retval <- rbind(bjones.ts(rets),wrap.itheta(rets,...))
	rownames(retval) <- c("Britten Jones t-stat","via itheta_vcov")
	return(retval)
}
do.all <- function(rets,...) {
	set.coln(rets)
	retval <- do.both(rets,...)
	retval <- rbind(retval,wrap.mp(rets,...))
	rownames(retval)[dim(retval)[1]] <- "via mp_vcov"
	return(retval)
}
n.day <- 1000
n.stock <- 5

## ----'me_vs_bjones',echo=TRUE---------------------------------
# under the null: all returns are zero mean;
set.seed(as.integer(charToRaw("7af85b0b-e521-4059-bebe-55ad9a9a0456")))
rets <- matrix(rnorm(n.day * n.stock),nrow=n.day)
# compare them:
print(do.all(rets))

# returns under the alternative
set.seed(as.integer(charToRaw("464dcbea-375b-49d0-8f38-b48b5a33c7ea")))
rets <- matrix(rnorm(n.day * n.stock,mean=1e-1),nrow=n.day)
print(do.all(rets))

## ----'get_mp',echo=TRUE---------------------------------------
# returns under the alternative
set.seed(as.integer(charToRaw("464dcbea-375b-49d0-8f38-b48b5a33c7ea")))
rets <- matrix(rnorm(n.day * n.stock,mean=1e-1),nrow=n.day)
print(wrap.itheta(rets,ret.what='mp'))

## ----'multiple_correlation',echo=TRUE,cache=FALSE-------------
# multiple correlation coefficients of portfolio
# error to precision errors.
mult.cor <- function(rets,...) {
	set.coln(rets)
	# 'flipping the sign on returns is idiomatic'
	asymv <- itheta_vcov(- as.matrix(rets),...)
	Ro <- cov2cor(asymv$Ohat)
	prec.idx <- (asymv$p+1):(dim(asymv$Ohat)[1])
	prec.Ro <- Ro[prec.idx,prec.idx]
	xcor <- Ro[2:asymv$p,prec.idx]
	R.sq <- diag(xcor %*% (solve(prec.Ro,t(xcor))))
}
set.seed(as.integer(charToRaw("464dcbea-375b-49d0-8f38-b48b5a33c7ea")))
rets <- matrix(rnorm(n.day * n.stock,mean=8e-2),nrow=n.day)
print(signif(100 * mult.cor(rets),digits=2))

## ----'multiple_correlation2',echo=TRUE,cache=FALSE------------
set.seed(as.integer(charToRaw("464dcbea-375b-49d0-8f38-b48b5a33c7ea")))
rets <- matrix(rnorm(n.day * n.stock,mean=1.6e-1),nrow=n.day)
print(signif(100 * mult.cor(rets),digits=2))

## ----'conditional',echo=TRUE----------------------------------
# generate data with given W, Sigma
Xgen <- function(W,Sigma,Feat) {
 Btrue <- Sigma %*% W
 Xmean <- Feat %*% t(Btrue)
 Shalf <- chol(Sigma)
 X <- Xmean + matrix(rnorm(prod(dim(Xmean))),ncol=dim(Xmean)[2]) %*% Shalf
}
n.feat <- 4
n.ret <- 8
n.obs <- 2000
set.seed(12321)

Feat <- matrix(rnorm(n.obs * n.feat),ncol=n.feat)
Wtrue <- 10 * matrix(rnorm(n.feat * n.ret),ncol=n.feat)
Sigma <- cov(matrix(rnorm(100*n.ret),ncol=n.ret))
Sigma <- Sigma + diag(seq(from=1,to=3,length.out=n.ret))
X <- Xgen(Wtrue,Sigma,Feat)
ism <- mp_vcov(X,feat=Feat,fit.intercept=TRUE)

# a bit of legerdemain b/c there's an intercept term
# fit
Wcomp <- cbind(0,Wtrue)

## ----'scatfit',echo=TRUE,cache=FALSE,fig.cap=paste("Scatter of the fit value against the true value of the Markowitz Coefficien is shown, for",n.ret,"assets, and",n.feat,"predictive features, given",n.obs,"days of observations of Gaussian returns. The $y=x$ line is plotted in green.")----
# scatter them against each other
w.true <- Wcomp
w.fit <- ism$W
dim(w.true) <- c(length(w.true),1)
dim(w.fit) <- c(length(w.fit),1)
plot(w.true, w.fit, main="Markowitz coefficient",
  	 xlab="True Value ", ylab="Fit Value ", pch=1)
#abline(lm(w.fit ~ w.true), col="red") 
abline(a=0,b=1, col="green") 

## ----'check_normality',echo=TRUE,cache=FALSE,fig.cap=paste("Sample quantiles of the error of the \\txtMC, transformed to approximate unit covariance using the estimated covariance, are plotted against those of the normal.")----
n.feat <- 4
n.ret <- 8
n.obs <- 2000
set.seed(12321)
Feat <- matrix(rnorm(n.obs * n.feat),ncol=n.feat)
Wtrue <- 10 * matrix(rnorm(n.feat * n.ret),ncol=n.feat)
Sigma <- cov(matrix(rnorm(100*n.ret),ncol=n.ret))
Sigma <- Sigma + diag(seq(from=1,to=3,length.out=n.ret))
X <- Xgen(Wtrue,Sigma,Feat)
ism <- mp_vcov(X,feat=Feat,fit.intercept=TRUE)
Wcomp <- cbind(0,Wtrue)
# compute the errors
errs <- ism$W - Wcomp
dim(errs) <- c(length(errs),1)
# transform them to approximately identity covariance
Zerr <- solve(t(chol(ism$What)),errs)
# did it work?
qqnorm(Zerr)
qqline(Zerr,col=2)

## ----'cons_zero',echo=TRUE,cache=FALSE------------------------
# first for the case where the real Markowitz Portfolio is actually
# just 'the market': equal dollar long in each stock.
set.seed(as.integer(charToRaw("dbeebe3f-da7e-4d11-b014-feac88a1d6cb")))
rets <- matrix(rnorm(n.day * n.stock,mean=1e-1),nrow=n.day)
Jmat <- matrix(c(1,1,1,-1),nrow=2,ncol=2*n.stock)  # overbuilt
Jmat <- Jmat[,1:n.stock]  # fixed.
print(Jmat)
# first, unconstrained:
print(wrap.mp(rets))
# now with subspace constraint:
print(wrap.mp(rets,Jmat=Jmat))
# and print the portfolio too:
print(wrap.mp(rets,ret.what='mp',Jmat=Jmat))

# now for the case where the real Markowitz portfolio is not just the market.
set.seed(as.integer(charToRaw("420f1dfd-b19b-4175-83b3-b96548264bf8")))
rets <- matrix(rnorm(n.day * n.stock,mean=0),nrow=n.day)
mu.off <- t(matrix(seq(from=-0.15,to=0.15,length.out=n.stock),nrow=n.stock,ncol=n.day))
rets <- rets + mu.off
print(wrap.mp(rets,Jmat=Jmat))
# and print the portfolio too:
print(wrap.mp(rets,ret.what='mp',Jmat=Jmat))

## ----'cons_one',echo=TRUE,cache=FALSE-------------------------

# first for the case where the real Markowitz Portfolio is actually
# equal dollar long in each stock.
set.seed(as.integer(charToRaw("0bda3ab6-53a7-4f5a-aa6a-0fe21edbaa20")))
rets <- matrix(rnorm(n.day * n.stock,mean=1e-1),nrow=n.day)
Gmat <- matrix(1,nrow=1,ncol=n.stock)
print(wrap.mp(rets,Gmat=Gmat))
# and print the portfolio too:
print(wrap.mp(rets,ret.what='mp',Gmat=Gmat))
# and in the case where it is not.
set.seed(as.integer(charToRaw("420f1dfd-b19b-4175-83b3-b96548264bf8")))
rets <- matrix(rnorm(n.day * n.stock,mean=0),nrow=n.day)
mu.off <- t(matrix(seq(from=-0.8,to=0.8,length.out=n.stock),nrow=n.stock,ncol=n.day))
rets <- rets + mu.off
Gmat <- matrix(1,nrow=1,ncol=n.stock)
print(wrap.mp(rets,Gmat=Gmat))
# and print the portfolio too:
print(wrap.mp(rets,ret.what='mp',Gmat=Gmat))

## ----'rao_giri',echo=TRUE,eval=TRUE---------------------------
rao.giri <- function(rets,Gmat,...) {
	set.coln(rets)
	asymv <- mp_vcov(as.matrix(rets),Gmat=Gmat,...)
	stat <- asymv$mu[1]
	a.var <- asymv$Ohat[1,1]
	return(list(stat=stat,a.var=a.var,wald=stat/sqrt(a.var)))
}
# here we hedge out G, the rows of which span the true
# Markowitz portfolio
set.seed(as.integer(charToRaw("4618fc2e-9c58-4ea7-81f7-6ed771ace500")))
rets <- matrix(rnorm(n.day * n.stock,mean=1e-1),nrow=n.day)
Gmat <- matrix(1,nrow=1,ncol=n.stock)
print(rao.giri(rets,Gmat)$wald)

# here we hedge out a random G, the rows of which do not
# span the true Markowitz portfolio
set.seed(as.integer(charToRaw("4abaccd9-8cac-4149-b30d-f3b4d32b44df")))
rets <- matrix(rnorm(n.day * n.stock,mean=1e-1),nrow=n.day)
Gmat <- matrix(rnorm(n.stock),nrow=1,ncol=n.stock)
print(rao.giri(rets,Gmat)$wald)

## ----'ff_load',eval=FALSE,echo=TRUE---------------------------
#  ff.data <- read.csv(paste0('http://www.quandl.com/api/v1/datasets/',
#  'KFRENCH/FACTORS_M.csv?&trim_start=1926-07-31&trim_end=2013-10-31',
#  '&sort_order=asc'), colClasses=c('Month'='Date'))
#  
#  rownames(ff.data) <- ff.data$Month
#  ff.data <- ff.data[,! (colnames(ff.data) %in% c("Month"))]
#  # will not matter, but convert pcts:
#  ff.data <- 1e-2 * ff.data

## ----'ff_load_sneaky',eval=TRUE,echo=FALSE--------------------
# sleight of hand to load precomputed data instead.
fname <- system.file('extdata','ff_data.rda',package='MarkowitzR')
if (fname == "") {
	fname <- 'ff_data.rda'
}
# poofs ff.data here
load(fname)

## ----'ff_process',echo=TRUE,cachc=TRUE------------------------
rfr <- ff.data[,'RF']

ff.ret <- cbind(ff.data[,'Mkt.RF'],
	ff.data[,c('HML','SMB')] - rep(rfr,2))
colnames(ff.ret)[1] <- "MKT"
# try both procedures:
print(do.both(ff.ret))
# try with robust s.e.
if (require(sandwich)) {
	print(do.both(ff.ret,vcov.func=sandwich::vcovHAC))
}

