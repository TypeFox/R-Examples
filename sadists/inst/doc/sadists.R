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
opts_chunk$set(cache=TRUE,cache.path="cache/sadists")

#opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
opts_chunk$set(fig.path="figure/sadists",dev=c("pdf"))
#opts_chunk$set(fig.width=4.0,fig.height=6,dpi=200)
# see http://stackoverflow.com/q/23419782/164611
opts_chunk$set(fig.width=5.0,fig.height=3.0,out.width="5in",out.height="3in",dpi=200)

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
#opts_chunk$set(size="tiny")
opts_chunk$set(size="small")
opts_chunk$set(tidy=TRUE,tidy.opts=list(width.cutoff=36,blank=TRUE))

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

library(sadists)

gen_norm <- rnorm
lseq <- function(from,to,length.out) { 
	exp(seq(log(from),log(to),length.out = length.out))
}

## ----'setup',echo=TRUE,eval=TRUE------------------------------
require(ggplot2)
require(grid)
testf <- function(dpqr,nobs,...) {
	rv <- sort(dpqr$r(nobs,...))
	data <- data.frame(draws=rv,pvals=dpqr$p(rv,...))
	text.size <- 6   # sigh

	# http://stackoverflow.com/a/5688125/164611
	p1 <- qplot(rv, geom = 'blank') +   
		geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +  
		stat_function(fun = function(x) { dpqr$d(x,...) }, aes(colour = 'Theoretical')) +                       
		geom_histogram(aes(y = ..density..), alpha = 0.3) +                        
		scale_colour_manual(name = 'Density', values = c('red', 'blue')) +
		theme(text=element_text(size=text.size)) + 
		labs(title="Density (tests dfunc)")

	# Q-Q plot
	p2 <- ggplot(data, aes(sample = draws)) + stat_qq(distribution=function(p) { dpqr$q(p,...) }) +
		geom_abline(slope=1,intercept=0,colour='red') + 
		theme(text=element_text(size=text.size)) + 
		labs(title="Q-Q plot (tests qfunc)")

	# empirical CDF of the p-values; should be uniform
	p3 <- ggplot(data, aes(sample = pvals)) + stat_qq(distribution=qunif) + 
		geom_abline(slope=1,intercept=0,colour='red') + 
		theme(text=element_text(size=text.size)) + 
		labs(title="P-P plot (tests pfunc)")

	# Define grid layout to locate plots and print each graph
	pushViewport(viewport(layout = grid.layout(2, 2)))
	print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
	print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
	print(p3, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
}

## ----'sumchisqpow',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the sum of chi-squares to a power distribution."----
require(sadists)
wts <- c(-1,1,3,-3)
df <- c(100,200,100,50)
ncp <- c(0,1,0.5,2)
pow <- c(1,0.5,2,1.5)
testf(list(d=dsumchisqpow,p=psumchisqpow,q=qsumchisqpow,r=rsumchisqpow),nobs=2^14,wts,df,ncp,pow)

## ----'kprime',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the K-prime distribution."----
require(sadists)
v1 <- 50
v2 <- 80
a <- 0.5
b <- 1.5
testf(list(d=dkprime,p=pkprime,q=qkprime,r=rkprime),nobs=2^14,v1,v2,a,b)

## ----'lambdap',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the Lambda-prime distribution."----
require(sadists)
df <- 50
ts <- 1.5
testf(list(d=dlambdap,p=plambdap,q=qlambdap,r=rlambdap),nobs=2^14,df,ts)

## ----'upsilon',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the upsilon distribution."----
require(sadists)
df <- c(30,50,100,20,10)
ts <- c(-3,2,5,-4,1)
testf(list(d=dupsilon,p=pupsilon,q=qupsilon,r=rupsilon),nobs=2^14,df,ts)

## ----'dnf',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the doubly non-central F distribution."----
require(sadists)
df1 <- 40
df2 <- 80
ncp1 <- 1.5
ncp2 <- 2.5
testf(list(d=ddnf,p=pdnf,q=qdnf,r=rdnf),nobs=2^14,df1,df2,ncp1,ncp2)

## ----'dnt',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the doubly non-central t distribution."----
require(sadists)
df <- 75
ncp1 <- 2
ncp2 <- 3
testf(list(d=ddnt,p=pdnt,q=qdnt,r=rdnt),nobs=2^14,df,ncp1,ncp2)

## ----'dnbeta',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the doubly non-central Beta distribution."----
require(sadists)
df1 <- 40
df2 <- 80
ncp1 <- 1.5
ncp2 <- 2.5
testf(list(d=ddnbeta,p=pdnbeta,q=qdnbeta,r=rdnbeta),nobs=2^14,df1,df2,ncp1,ncp2)

## ----'dneta',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the doubly non-central Eta distribution."----
require(sadists)
df <- 100
ncp1 <- 0.5
ncp2 <- 2.5
testf(list(d=ddneta,p=pdneta,q=qdneta,r=rdneta),nobs=2^14,df,ncp1,ncp2)

## ----'sumlogchisq',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the sum of log of chi-squares distribution."----
require(sadists)
wts <- c(5,-4,10,-15)
df <- c(100,200,100,50)
ncp <- c(0,1,0.5,2)
testf(list(d=dsumlogchisq,p=psumlogchisq,q=qsumlogchisq,r=rsumlogchisq),nobs=2^14,wts,df,ncp)

## ----'prodchisqpow',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the product of chi-squares to a power distribution."----
require(sadists)
df <- c(100,200,100,50)
ncp <- c(0,1,0.5,2)
pow <- c(1,0.5,2,1.5)
testf(list(d=dprodchisqpow,p=pprodchisqpow,q=qprodchisqpow,r=rprodchisqpow),nobs=2^14,df,ncp,pow)

## ----'proddnf',echo=TRUE,eval=TRUE,fig.cap="Confirming the dpqr functions of the product of doubly non-central Fs distribution."----
require(sadists)
df1 <- c(10,20,5)
df2 <- c(1000,500,150)
ncp1 <- c(1,0,2.5)
ncp2 <- c(0,1.5,5)
testf(list(d=dproddnf,p=pproddnf,q=qproddnf,r=rproddnf),nobs=2^14,df1,df2,ncp1,ncp2)

