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
opts_chunk$set(cache=TRUE,cache.path="cache/asymarko")

#opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
opts_chunk$set(fig.path="figure/asymarko",dev=c("pdf"))
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

# library(quantmod)
# options("getSymbols.warning4.0"=FALSE)

library(MarkowitzR)

gen_norm <- rnorm
lseq <- function(from,to,length.out) { 
	exp(seq(log(from),log(to),length.out = length.out))
}

