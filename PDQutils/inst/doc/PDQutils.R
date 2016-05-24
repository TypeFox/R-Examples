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
opts_chunk$set(cache=TRUE,cache.path="cache/PDQutils")

#opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
opts_chunk$set(fig.path="figure/PDQutils",dev=c("pdf"))
opts_chunk$set(fig.width="4in",fig.height="4in",fig.width=8,fig.height=8,dpi=450)

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
opts_chunk$set(tidy=TRUE,tidy.opts=list(width.cutoff=50,blank=TRUE))

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

library(PDQutils)

gen_norm <- rnorm
lseq <- function(from,to,length.out) { 
	exp(seq(log(from),log(to),length.out = length.out))
}

## ----'rsnak',eval=TRUE,echo=TRUE------------------------------
rsnak <- function(n,dfs) {
	samples <- Reduce('+',lapply(dfs,function(k) { sqrt(rchisq(n,df=k)/k) }))
}

## ----'testit',eval=TRUE,echo=TRUE,fig.width=4.5,fig.height=4,fig.scap=NA,fig.cap=paste0("A q-q plot of ",n.samp," draws from a sum of Nakagamis distribution is shown against normality.")----
n.samp <- 1e5
dfs <- c(8,15,4000,10000)
set.seed(18181)
# now draw from the distribution 
rvs <- rsnak(n.samp,dfs)
data <- data.frame(draws=rvs)

library(ggplot2)
mu <- mean(rvs)
sigma <- sd(rvs)
ph <- ggplot(data, aes(sample = draws)) + stat_qq(distribution=function(p) { qnorm(p,mean=mu,sd=sigma) }) +
		geom_abline(slope=1,intercept=0,colour='red') + 
		theme(text=element_text(size=8)) + 
		labs(title="Q-Q plot (against normality)")

print(ph)

## ----'snakcu',eval=TRUE,echo=TRUE-----------------------------
# for the moment2cumulant function:
library(PDQutils)

# compute the first ord.max raw cumulants of the sum of Nakagami variates
snak_cumulants <- function(dfs,ord.max=10) {
	# first compute the raw moments
	moms <- lapply(dfs,function(nu) { 
		ords <- 1:ord.max
		moms <- 2 ^ (ords/2.0) * exp(lgamma((nu+ords)/2) - lgamma(nu/2))
		# we are dividing the chi by sqrt the d.f.
		moms <- moms / (nu ^ (ords/2.0))
		moms
	})
	# turn moments into cumulants
	cumuls <- lapply(moms,moment2cumulant)
	# sum the cumulants
	tot.cumul <- Reduce('+',cumuls)
	return(tot.cumul)
}

## ----'dpqsnak',eval=TRUE,echo=TRUE----------------------------

library(PDQutils)

dsnak <- function(x,dfs,ord.max=10,...) {
	raw.moment <- cumulant2moment(snak_cumulants(dfs,ord.max))
	retval <- dapx_gca(x,raw.moment,support=c(0,Inf),...)
	return(retval)
}
psnak <- function(q,dfs,ord.max=10,...) {
	raw.moment <- cumulant2moment(snak_cumulants(dfs,ord.max))
	retval <- papx_gca(q,raw.moment,support=c(0,Inf),...)
	return(retval)
}
qsnak <- function(p,dfs,ord.max=10,...) {
	raw.cumul <- snak_cumulants(dfs,ord.max)
	retval <- qapx_cf(p,raw.cumul,support=c(0,Inf),...)
	return(retval)
}

## ----'dpqsnak2',eval=TRUE,echo=TRUE---------------------------

dsnak_2 <- function(x,dfs,ord.max=10,...) {
	raw.cumul <- snak_cumulants(dfs,ord.max)
	retval <- dapx_edgeworth(x,raw.cumul,support=c(0,Inf),...)
	return(retval)
}
psnak_2 <- function(q,dfs,ord.max=10,...) {
	raw.cumul <- snak_cumulants(dfs,ord.max)
	retval <- papx_edgeworth(q,raw.cumul,support=c(0,Inf),...)
	return(retval)
}

## ----'improvedqq',eval=TRUE,echo=TRUE,fig.width=4.5,fig.height=4,fig.scap=NA,fig.cap=paste0("A q-q plot of ",n.samp," draws from a sum of Nakagamis distribution is shown against quantiles from the `qsnak' function.")----

data <- data.frame(draws=rvs)
library(ggplot2)
ph <- ggplot(data, aes(sample = draws)) + stat_qq(distribution=function(p) { qsnak(p,dfs=dfs) }) +
		geom_abline(slope=1,intercept=0,colour='red') + 
		theme(text=element_text(size=8)) + 
		labs(title="Q-Q against qsnak (C-F appx.)")
print(ph)


## ----'psnak_ecdf',eval=TRUE,echo=TRUE,fig.width=4.5,fig.height=4,fig.scap=NA,fig.cap=paste0("The empirical CDF of the approximate CDF of a sum of Nakagamis distribution on ",n.samp," draws is shown.")----

apx.p <- psnak(rvs,dfs=dfs)
require(ggplot2)
ph <- ggplot(data.frame(pv=apx.p),aes(x=pv)) + stat_ecdf(geom='step')
print(ph)

## ----'chisetup',eval=TRUE,echo=TRUE,fig.width=4.5,fig.height=4,fig.scap=NA,fig.cap="The true PDF of a normalized $\\chi^2_5$ distribution is shown, along with the 2- and 6-term Gram Charlier approximations. This replicates Figure 1 of Blinnikov and Moessner. \\cite{1998A&AS..130..193B}"----
# compute moments and cumulants:
df <- 5
max.ord <- 20
subords <- 0:(max.ord - 1)
raw.cumulants <- df * (2^subords) * factorial(subords)
raw.moments <- cumulant2moment(raw.cumulants)

# compute the PDF of the rescaled variable:
xvals <- seq(-sqrt(df/2) * 0.99,6,length.out=1001)
chivals <- sqrt(2*df) * xvals + df
pdf.true <- sqrt(2*df) * dchisq(chivals,df=df)

pdf.gca2 <- sqrt(2*df) * dapx_gca(chivals,raw.moments=raw.moments[1:2],support=c(0,Inf))
pdf.gca6 <- sqrt(2*df) * dapx_gca(chivals,raw.moments=raw.moments[1:6],support=c(0,Inf))

all.pdf <- data.frame(x=xvals,true=pdf.true,gca2=pdf.gca2,gca6=pdf.gca6)

# plot it by reshaping and ggplot'ing
require(reshape2)
arr.data <- melt(all.pdf,id.vars='x',variable.name='pdf',value.name='density')

require(ggplot2)
ph <- ggplot(arr.data,aes(x=x,y=density,group=pdf,colour=pdf)) + geom_line()
print(ph)

## ----'chitwo',eval=TRUE,echo=TRUE,fig.width=4.5,fig.height=4,fig.scap=NA,fig.cap="The true PDF of a normalized $\\chi^2_5$ distribution is shown, along with the 2- and 4-term Edgeworth expansions. This replicates Figure 6 of Blinnikov and Moessner. \\cite{1998A&AS..130..193B}"----
# compute the PDF of the rescaled variable:
xvals <- seq(-sqrt(df/2) * 0.99,6,length.out=1001)
chivals <- sqrt(2*df) * xvals + df
pdf.true <- sqrt(2*df) * dchisq(chivals,df=df)

pdf.edgeworth2 <- sqrt(2*df) * dapx_edgeworth(chivals,raw.cumulants=raw.cumulants[1:4],support=c(0,Inf))
pdf.edgeworth4 <- sqrt(2*df) * dapx_edgeworth(chivals,raw.cumulants=raw.cumulants[1:6],support=c(0,Inf))

all.pdf <- data.frame(x=xvals,true=pdf.true,edgeworth2=pdf.edgeworth2,edgeworth4=pdf.edgeworth4)

# plot it by reshaping and ggplot'ing
require(reshape2)
arr.data <- melt(all.pdf,id.vars='x',variable.name='pdf',value.name='density')

require(ggplot2)
ph <- ggplot(arr.data,aes(x=x,y=density,group=pdf,colour=pdf)) + geom_line()
print(ph)

