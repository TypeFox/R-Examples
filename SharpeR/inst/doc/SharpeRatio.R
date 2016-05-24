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
opts_chunk$set(cache=TRUE,cache.path="cache/SharpeRatio")

#opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
opts_chunk$set(fig.path="figure/SharpeRatio",dev=c("pdf"))
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

mc.resolution <- ifelse(LONG.FORM,1000,200)
mc.resolution <- max(mc.resolution,100)

library(quantmod)
options("getSymbols.warning4.0"=FALSE)

library(SharpeR)

gen_norm <- rnorm
lseq <- function(from,to,length.out) { 
	exp(seq(log(from),log(to),length.out = length.out))
}

## ----'generate_tbias',include=FALSE---------------------------
# the bias of the t
f_tbias <- function(n) { 
	sqrt((n-1) / 2) * exp(lgamma((n-2)/2) - lgamma((n-1)/2))
}
# approximate tbias
f_apx_tbias1 <- function(n) { 
	1 + (0.75 / n)
}
f_apx_tbias2 <- function(n) { 
	#1 + (0.75 / n) + (2 / n**2)
	1 + (0.75 / (n - 1)) + (32 / (25 * ((n-1) ** 2)))
}
n <- 12
tbias <- f_tbias(n)

## ----'f_tpower',echo=FALSE------------------------------------
# the power of the univariate t-test;
f_tpower <- function(n,rho = 0,alpha = 0.05) {
	f_tpower_ncp(ncp = sqrt(n) * rho,n = n,alpha = alpha)
}

# the power of the univariate t-test as function of the ncp
f_tpower_ncp <- function(ncp,n,alpha = 0.05) {
	pt(qt(1-alpha,n-1),df = n-1,ncp = ncp,lower.tail = FALSE)
}

# the power of an f-test 
f_fpower <- function(df1,df2,ncp,alpha = 0.05) {
	pf(qf(alpha,df1=df1,df2=df2,ncp=0,lower.tail=FALSE),df1 = df1,df2 = df2,ncp = ncp,lower.tail = FALSE)
}

# find the sample size for a given power of the univariate t-test
f_treqsize <- function(rho,powr = 0.80,alpha = 0.05) {
	zz <- uniroot(function(n,rho,alpha,powr)(f_tpower(n,rho,alpha) - powr),
								c(8,20 * 10 / (rho*rho)),rho = rho,powr = powr,alpha = alpha)
	return(zz$root)
}

## ----'sr_power_table',echo=FALSE,results='asis'---------------
dpy <- 253
rhovals <- seq(0.5,3.0,by=0.10)
rhovals_day <- rhovals / sqrt(dpy)
samps25.1 <- (1/dpy) * unlist(lapply(rhovals_day, f_treqsize, powr = 0.25))
samps50.1 <- (1/dpy) * unlist(lapply(rhovals_day, f_treqsize, powr = 0.50))
samps80.1 <- (1/dpy) * unlist(lapply(rhovals_day, f_treqsize, powr = 0.80))
#2sided test
samps25.2 <- (1/dpy) * unlist(lapply(rhovals_day, f_treqsize, alpha=0.025, powr = 0.25))
samps50.2 <- (1/dpy) * unlist(lapply(rhovals_day, f_treqsize, alpha=0.025, powr = 0.50))
samps80.2 <- (1/dpy) * unlist(lapply(rhovals_day, f_treqsize, alpha=0.025, powr = 0.80))
foo <- data.frame("one sided" = c(median(samps25.1 * rhovals**2),median(samps50.1 * rhovals**2),median(samps80.1 * rhovals**2)),
									"two sided" = c(median(samps25.2 * rhovals**2),median(samps50.2 * rhovals**2),median(samps80.2 * rhovals**2)),
									row.names = c("power = 0.25","power = 0.50","power = 0.80"))
library(xtable)
xres <- xtable(foo,label="tab:ttestpower",caption="Scaling of sample size with respect to $\\psnr^2$ required to achieve a fixed power in the t-test, at a fixed $\\typeI = 0.05$ rate.")
print(xres,include.rownames=TRUE)

## ----'power_thing',echo=FALSE,fig.cap="The percent error of the power mnemonic $e\\approx\\ssiz \\psnrsq$ is plotted versus \\psnr."----
ope <- 253
 zetas <- seq(0.1,2.5,length.out=51)
ssizes <- sapply(zetas,function(zed) { 
 x <- power.sr_test(n=NULL,zeta=zed,sig.level=0.05,power=0.5,ope=ope)
 x$n / ope 
})
plot(zetas,100 * ((exp(1) / zetas^2) - ssizes)/ssizes, ylab="error in mnemonic rule (as %)")

## ----'sobering',include=FALSE---------------------------------
foo.power <- power.sr_test(n=253,zeta=NULL,sig.level=0.05,power=0.5,ope=253)

## ----"autocorr_setup",echo=FALSE------------------------------
fname <- system.file('extdata','autocorr_study.rda',package='SharpeR')
if (fname == "") {
	fname <- 'autocorr_study.rda'
}
# poofs phivals empiricals predicteds ntrials nyr
load(fname)

# maybe need these some other time?
# vanilla Sharpe ratio in terms of whatever input units
f_vsharpe <- function(rets) {
	return(mean(rets) / sd(rets))
}
# annualized Sharpe ratio
f_asharpe <- function(rets, dpy = 253) {
	return(sqrt(dpy) * f_vsharpe(rets))
}
# the t statistic
f_tstat <- function(rets) { 
	return(sqrt(length(rets)) * f_vsharpe(rets))
}

## ----'autocorr_spread',echo=FALSE,fig.cap=paste("The empirical standard deviation for the \\tstat-statistic is shown at different values of the autocorrelation, \\pacor. Each point represents",ntrials,"series of approximately",nyr,"years of daily data, with each series generated by an AR(1) process with normal innovations.  Each series has actual SNR of zero. The fit line is that suggested by Van Belle's correction for autocorrelation, namely $\\sqrt{(1 + \\pacor) / (1 - \\pacor)}$.")----
plot(phivals,empiricals,pch=1,col="red",
	ylab="standard deviation of t statistic",
	xlab="autocorrelation")
lines(phivals,predicteds,col="black")
legend("topleft",c("Empirical","Predicted"),
	lty=c(0,1),pch=c(1,NA_integer_),col=c("red","black"))

## ----'leverage_foo',include=FALSE-----------------------------
Elevi <- 1.5 
Elevis <- Elevi ** 2
Vlevi <- 1 / 12
ratlevi <- Vlevi / Elevis
vixlevi <- (0.4)^2

## ----'skewstudy_prelim',echo=FALSE,cache=FALSE----------------
# MOCK it up.
fname <- system.file('extdata','skew_study.rda',package='SharpeR')
if (fname == "") {
	fname <- 'skew_study.rda'
}
# poofs res, ntrials, SPX.rets, dpy
load(fname)

SPX.start <- time(SPX.rets[1])
SPX.end <- time(SPX.rets[length(SPX.rets)])

## ----'skewstudy',echo=FALSE,results='asis',cache=FALSE--------
# res, ntrials are poofed above.
#digits=c(0,0,0,2,2,2),
xres <- xtable(res,label="tab:sharpe_skew_robustness",
							 display=c('s','s','s','g','g','fg','fg'),
							 align="ll|lrr|rr",
							 caption=paste("Empirical type I rates of the test for $\\psnr = 1.0\\yrtomhalf$ via distribution of the \\txtSR are given for various distributions of returns.  The empirical rates are based on ", 
														 ntrials, 
														 " simulations of three years of daily returns, with a nominal rate of $\\typeI = 0.05$. The 'corrected' type I rates refer to a normal approximation using Mertens' correction. Skew appears to have a much more adverse effect than kurtosis alone."))
print(xres,include.rownames=FALSE,hline.after=c(0,0,2,4,7,10,dim(res)[1]))

## ----'hot_ssiz_table',echo=FALSE,results='asis'---------------
# MOCK it up.
fname <- system.file('extdata','hotelling_power_rule.rda',package='SharpeR')
if (fname == "") {
	fname <- 'hotelling_power_rule.rda'
}
# poofs res, ntrials, SPX.rets, dpy
load(fname)
xres <- xtable(hot.power,label="tab:htestpower",
	caption="The numerator in the sample size relationship required to achieve a fixed power in Hotelling's test is shown. The type I rate is 0.05.")
print(xres,include.rownames=TRUE,sanitize.text.function=function(x){x})

## ----'example_usage_Fpower',include=FALSE---------------------
ex_p <- 30
ex_n <- 1000
ex_snr <- 1.5
dpy <- 253
ex_snr_d <- ex_snr / sqrt(dpy)
c <- ex_p / ex_n
cplus <- (c + ex_snr_d ** 2)
meanv <- sqrt(cplus / (1 - c))
stdv <- sqrt((1 / ex_n) * (ex_snr_d ** 4 + 2 * ex_snr_d ** 2 + c) / (2 * cplus * (1 - c) ** 2))

## ----'haircut_study_load',echo=FALSE,cache=FALSE--------------
# MOCK it up.
fname <- system.file('extdata','haircut_study.rda',package='SharpeR')
if (fname == "") {
	fname <- 'haircut_study.rda'
}
# poofs n.sim,n.stok,n.yr,n.obs,zeta.s,ope,hcuts
load(fname)

medv.true <- median(hcuts)
med.snr.true <- zeta.s * (1 - medv.true)

## ----'haircut_study_mock',eval=FALSE,echo=TRUE----------------
#  require(MASS)
#  
#  # simple markowitz.
#  simple.marko <- function(rets) {
#  	mu.hat <- as.vector(apply(rets,MARGIN=2,mean,na.rm=TRUE))
#  	Sig.hat <- cov(rets)
#  	w.opt <- solve(Sig.hat,mu.hat)
#  	retval <- list('mu'=mu.hat,'sig'=Sig.hat,'w'=w.opt)
#  	return(retval)
#  }
#  # make multivariate pop. & sample w/ given zeta.star
#  gen.pop <- function(n,p,zeta.s=0) {
#  	true.mu <- matrix(rnorm(p),ncol=p)
#  	#generate an SPD population covariance. a hack.
#  	xser <- matrix(rnorm(p*(p + 100)),ncol=p)
#  	true.Sig <- t(xser) %*% xser
#  	pre.sr <- sqrt(true.mu %*% solve(true.Sig,t(true.mu)))
#  	#scale down the sample mean to match the zeta.s
#  	true.mu <- (zeta.s/pre.sr[1]) * true.mu
#    X <- mvrnorm(n=n,mu=true.mu,Sigma=true.Sig)
#  	retval = list('X'=X,'mu'=true.mu,'sig'=true.Sig,'SNR'=zeta.s)
#  	return(retval)
#  }
#  # a single simulation
#  sample.haircut <- function(n,p,...) {
#  	popX <- gen.pop(n,p,...)
#  	smeas <- simple.marko(popX$X)
#  	# I have got to figure out how to deal with vectors...
#  	ssnr <- (t(smeas$w) %*% t(popX$mu)) / sqrt(t(smeas$w) %*% popX$sig %*% smeas$w)
#  	hcut <- 1 - (ssnr / popX$SNR)
#  	# for plugin estimator, estimate zeta.star
#  	asro <- sropt(z.s=sqrt(t(smeas$w) %*% smeas$mu),df1=p,df2=n)
#  	zeta.hat.s <- inference(asro,type="KRS")  # or 'MLE', 'unbiased'
#  	return(c(hcut,zeta.hat.s))
#  }
#  
#  # set everything up
#  set.seed(as.integer(charToRaw("496509a9-dd90-4347-aee2-1de6d3635724")))
#  ope <- 253
#  n.sim <- 4096
#  n.stok <- 6
#  n.yr <- 4
#  n.obs <- ceiling(ope * n.yr)
#  zeta.s <- 1.20 / sqrt(ope)   # optimal SNR, in daily units
#  
#  # run some experiments
#  experiments <- replicate(n.sim,sample.haircut(n.obs,n.stok,zeta.s))
#  hcuts <- experiments[1,]

## ----'haircutting',fig.cap=paste("Q-Q plot of",n.sim,"simulated haircut values versus the approximation given by \\eqnref{hcut_apx} is shown.")----
print(summary(hcuts))
# haircut approximation in the equation above
qhcut <- function(p, df1, df2, zeta.s, lower.tail=TRUE) {
	atant <- atan((1/sqrt(df1-1)) * 
		qt(p,df=df1-1,ncp=sqrt(df2)*zeta.s,lower.tail=!lower.tail))
	# a slightly better approximation is:
	# retval <- 1 - sin(atant - 0.0184 * zeta.s * sqrt(df1 - 1))
	retval <- 1 - sin(atant)
}
# if you wanted to look at how bad the plug-in estimator is, then
# uncomment the following (you are warned):
# zeta.hat.s <- experiments[2,];                                   
# qqplot(qhcut(ppoints(length(hcuts)),n.stok,n.obs,zeta.hat.s),hcuts,
# 			 xlab = "Theoretical Approximate Quantiles", ylab = "Sample Quantiles");
# qqline(hcuts,datax=FALSE,distribution = function(p) { qhcut(p,n.stok,n.obs,zeta.hat.s) },
# 			 col=2)

# qqplot;
qqplot(qhcut(ppoints(length(hcuts)),n.stok,n.obs,zeta.s),hcuts,
			 xlab = "Theoretical Approximate Quantiles", ylab = "Sample Quantiles")
qqline(hcuts,datax=FALSE,distribution = function(p) { qhcut(p,n.stok,n.obs,zeta.s) },
			 col=2)

## ----'hcut_moments',echo=FALSE,results='asis'-----------------
mc.med <- median(hcuts)
mc.mean <- mean(hcuts)
mc.sd <- sd(hcuts)

mc.p <- n.stok
mc.n <- n.obs
mc.zs <- zeta.s

fit.med <- 1 - sin(atan(mc.zs * sqrt(mc.n / (mc.p - 1))))
fit.mean <- 1 - sqrt(1 - (mc.p / (mc.p + mc.n * mc.zs^2)))
fit.sd <- sqrt(mc.p) / (mc.p + (mc.n * mc.zs^2) ^ 1.08)

fit.df <- data.frame(Monte.Carlo=c(mc.med,mc.mean,mc.sd),
		approximation=c(fit.med,fit.mean,fit.sd))
rownames(fit.df) <- c("median","mean","standard deviation")

# 2FIX: start here;yy
xres <- xtable(fit.df,label="tab:hcutfit",
	caption=paste("Empirical approximate values of the median, mean, and",
	"standard deviation of the haircut distribution are given for",
	n.sim,"Monte Carlo simulations of",n.obs,"days of Gaussian data for",n.stok,
	"assets with $\\psnropt=",zeta.s * sqrt(ope),"\\yrtomhalf$.",
	"The approximations from \\eqnref{hcut_moment_appx} are also reported."))

print(xres,include.rownames=TRUE,sanitize.text.function=function(x){x})



## ----'pos_cons_example', include=FALSE, warning=FALSE, message=FALSE----
base.opt.1 <- function(rho,c1,c2) {
# compute v1, v2 that solve
#
#   min      v1^2 + v2^2
#   st.         v1 >= c1
#     -rho v1 + v2 >= c2 
# 
# return as a list?
# the prospective solutions are
# (0,0), (c1,0), (-2rho c2/(1+rho^2),c2(1-rho^2)/(1+rho^2)),
# and (c1,c2+rho*c1)
	rho2 <- rho^2
	oneprho2 <- 1 + rho2;
	solns <- matrix(c(0,0, c1,0, -2*rho*c2/(oneprho2),c2*(1-rho^2)/(oneprho2),
c1,c2+rho*c1),nrow=2)
	feasible <- solns[1,] >= c1
	solns <- solns[,feasible]
	feasible <- -rho * solns[1,] + solns[2,] >= c2
	solns <- solns[,feasible]
	solns <- matrix(solns,nrow=2)

	vals <- colSums(solns ^ 2)
	retval <- solns[,which.min(vals)]
	return(retval)
}
base.opt.2 <- function(rho,psi1,psi2) {
# compute lam1, lam2 that solve
# 
#   min     [lam1,lam2] * R^-1 [lam1,lam2]'
#   st.     R^-1 ([psi1,psi2]' + [lam1,lam2]') >= 0
#
#   where   R = [1,rho]
#               [rho,1]
#    so  R^-1 = [1, -rho]
#               [-rho, 1] / (1-rho^2)
#
# and then return R^-1 ([psi1,psi2]' + [lam1,lam2]')
	c1 <- rho * psi2 - psi1;
	c2 <- rho * psi1 - psi2;
	rv <- base.opt.1(rho,c1,c2)
	lam2 <- rv[2] / (1-rho^2)
	lam1 <- rv[1] + rho * lam2
	nmu1 <- psi1 + lam1
	nmu2 <- psi2 + lam2
	# deal with roundoff
	w1 <- max(0,nmu1 - rho * nmu2)
	w2 <- max(0,-rho * nmu1 + nmu2)
	return(c(w1,w2))
}
base.opt.3 <- function(mu1,mu2,sig1,sig12,sig2) {
# compute w1, w2 to solve
#
#   max   [w1,w2] * [mu1,mu2]'  / sqrt([w1,w2] * Sig * [w1,w2]')
#   st.   w1 >= 0
#         w2 >= 0
#
#  where  Sig = [sig1   sig12]
#               [sig12   sig2]
	rho <- sig12 / sqrt(sig1 * sig2)
	psi1 <- mu1 / sqrt(sig1)
	psi2 <- mu2 / sqrt(sig2)
	rv <- base.opt.2(rho,psi1,psi2)
	return(rv)
}
# test it
base.opt.3(1,0.5,1,0.4,1);
base.opt.3(1,0.5,1,0.5,1);
base.opt.3(1,0.5,1,0.6,1);
base.opt.3(1,0.5,1,0.8,1);

base.opt.3(1,0.5,1,-0.8,1);
opt.pos.T2 <- function(mu1,mu2,sig1,sig12,sig2) {
# compute the maximum value of 
#
#   max   [w1,w2] * [mu1,mu2]'  / sqrt([w1,w2] * Sig * [w1,w2]')
#   st.   w1 >= 0
#         w2 >= 0
	rv <- base.opt.3(mu1,mu2,sig1,sig12,sig2)
	w <- as.vector(rv)
	mu <- as.vector(c(mu1,mu2))
	Sig <- matrix(c(sig1,sig12,sig12,sig2),nrow=2)
	T <- (mu %*% w) / sqrt(t(w) %*% Sig %*% w)
	return(T^2)
}

	



## ----'me_vs_bjones',echo=TRUE---------------------------------
nday <- 1024
nstk <- 5

# under the null: all returns are zero mean;
set.seed(as.integer(charToRaw("7fbb2a84-aa4c-4977-8301-539e48355a35")))
rets <- matrix(rnorm(nday * nstk),nrow=nday)

# t-stat via Britten-Jones procedure
bjones.ts <- function(rets) {
	ones.vec <- matrix(1,nrow=dim(rets)[1],ncol=1)
	bjones.mod <- lm(ones.vec ~ rets - 1)
	bjones.sum <- summary(bjones.mod)
	retval <- bjones.sum$coefficients[,3]
}
# wald stat via inverse second moment trick
ism.ws <- function(rets) {
# flipping the sign on returns is idiomatic,
	asymv <- ism_vcov(- rets)
	asymv.mu <- asymv$mu[1:asymv$p]
	asymv.Sg <- asymv$Ohat[1:asymv$p,1:asymv$p]
	retval <- asymv.mu / sqrt(diag(asymv.Sg))
}

bjones.tstat <- bjones.ts(rets)
ism.wald <- ism.ws(rets)

# compare them:
print(bjones.tstat)
print(ism.wald)

# repeat under the alternative;
set.seed(as.integer(charToRaw("a5f17b28-436b-4d01-a883-85b3e5b7c218")))
zero.rets <- t(matrix(rnorm(nday * nstk),nrow=nday))
mu.vals <- (1/sqrt(253)) * seq(-1,1,length.out=nstk) 
rets <- t(zero.rets + mu.vals)

bjones.tstat <- bjones.ts(rets)
ism.wald <- ism.ws(rets)

# compare them:
print(bjones.tstat)
print(ism.wald)

