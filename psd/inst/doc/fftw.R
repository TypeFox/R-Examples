## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  rm(list=ls())
#  library(fftw)
#  library(rbenchmark)
#  library(plyr)
#  library(reshape2)
#  library(ggplot2)

## ----eval=FALSE, echo=TRUE, label="Benchmark function"-------------------
#  reps <- 10
#  dftbm <- function(nd, repls=reps){
#  	set.seed(1234)
#  	x <- rnorm(nd, mean=0, sd=1)
#  	bmd <- benchmark(replications=repls, fftw::FFT(x), stats::fft(x))
#  	bmd$num_dat <- nd
#  	bmd$relative[is.na(bmd$relative)] <- 1   # NA happens.
#  	return(bmd)
#  }

## ----eval=TRUE, echo=TRUE, label="Num. terms to bench."------------------
(nterms.even <- round(2**seq.int(from=4,to=20,by=1)))

## ----eval=FALSE, echo=TRUE, label="Use lapply to do benching."-----------
#  bench.even <- function(){
#    benchdat.e <- plyr::ldply(lapply(X=nterms.even, FUN=dftbm))
#    }
#  bench.even()

## ----eval=FALSE, echo=TRUE, label="Setup non highly composite lengths."----
#  nterms.odd <- nterms.even + 1
#  nterms.odd <- nterms.odd[nterms.odd < 50e3] # painfully long otherwise!

## ----eval=FALSE, echo=TRUE, label="Do benching."-------------------------
#  bench.odd <- function(){
#    benchdat.o <- plyr::ldply(lapply(X=nterms.odd, FUN=dftbm))
#    }
#  bench.odd() # FAIR WARNING: this can take a while!!

## ----eval=FALSE, echo=TRUE, label="Map/Reduce/Summarize etc"-------------
#  pltbench <- function(lentyp=c("even","odd")){
#    benchdat <- switch(match.arg(lentyp), even=benchdat.e, odd=benchdat.o)
#    stopifnot(exists("benchdat"))
#    tests <- unique(benchdat$test)
#    ## subset only information we care about
#    allbench.df.drp <- subset(benchdat,
#          select=c(test, num_dat, user.self, sys.self, elapsed, relative))
#    ## reduce data.frame with melt
#    allbench.df.mlt <- reshape2::melt(allbench.df.drp,
#                                      id.vars=c("test","num_dat"))
#    ## calculate the summary information to be plotted:
#    tmpd <- plyr::ddply(allbench.df.mlt,
#                        .(variable,  num_dat),
#                        summarise,
#                        summary="medians",
#                        value=ggplot2::mean_cl_normal(value)[1,1])
#    ## create copies for each test and map to data.frame
#    allmeds <<- plyr::ldply(lapply(X=tests,
#                                   FUN=function(x,df=tmpd){
#                                         df$test <- x; return(df)
#                                       }))
#    ## plot the benchmark data
#    # 1/sqrt(n) standard errors [assumes N(0,1)]
#    g <- ggplot(data=allbench.df.mlt,
#                aes(x=log10(num_dat),
#                    y=log2(value),
#                    ymin=log2(value*(1-1/sqrt(reps))),
#                    ymax=log2(value*(1+1/sqrt(reps))),
#                    colour=test,
#                    group=test)) +
#         scale_colour_discrete(guide="none") +
#         theme_bw()+
#         ggtitle(sprintf("DFT benchmarks of %s length series",toupper(lentyp))) +
#         ylim(c(-11,11))+
#         xlim(c(0.5,6.5))
#  
#    ## add previous summary curves if exist
#    if (exists("allmeds.prev")){
#       g <- g + geom_path(size=1.5, colour="dark grey", data=allmeds.prev,
#                          aes(group=test))
#                          }
#    ## create a facetted version
#    g2 <- g + facet_grid(variable~test) #, scales="free_y")
#    ## add the summary data as a line
#    g3 <- g2 + geom_path(colour="black", data=allmeds, aes(group=test))
#    ## and finally the data
#    print(g4 <<- g3 + geom_pointrange())
#  }

## ----eval=FALSE, echo=TRUE,  label="Plot non highly composite results."----
#  pltbench("even")
#  allmeds.prev <- allmeds
#  pltbench("odd")

## ----eval=TRUE, echo=TRUE, label=SI--------------------------------------
utils::sessionInfo()

