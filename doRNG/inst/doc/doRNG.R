
## ----pkgmaker_preamble, echo=FALSE, results='asis'-----------------------
pkgmaker::latex_preamble()


## ----bibliofile, echo=FALSE, results='asis'------------------------------
pkgmaker::latex_bibliography('doRNG')


## ----init, include = FALSE-----------------------------------------------
options(width=90)
library(pkgmaker)
library(knitr)
opts_chunk$set(size = "footnotesize")
knit_hooks$set(try = pkgmaker::hook_try)


## ----foreach-------------------------------------------------------------
# load and register parallel backend for multicore computations
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)

# perform 5 tasks in parallel
x <- foreach(i=1:5) %dopar% { 
	i + runif(1) 
}
unlist(x)


## ----dopar, tidy=FALSE---------------------------------------------------
# with standard %dopar%: foreach loops are not reproducible
set.seed(123)
res <- foreach(i=1:5) %dopar% { runif(3) }
set.seed(123)
res2 <- foreach(i=1:5) %dopar% { runif(3) }
identical(res, res2)


## ----dorng, tidy=FALSE---------------------------------------------------

# load the doRNG package
library(doRNG)

# using %dorng%: loops _are_ reproducible
set.seed(123)
res <- foreach(i=1:5) %dorng% { runif(3) }
set.seed(123)
res2 <- foreach(i=1:5) %dorng% { runif(3) }
identical(res, res2)


## ----attr----------------------------------------------------------------
attr(res, 'rng')


## ----nextRNGstream-------------------------------------------------------
nextRNGStream(c(407L, 1:6))


## ----options_single------------------------------------------------------
# use a single numeric as a seed
s <- foreach(i=1:5, .options.RNG=123) %dorng% { runif(3) }
s2 <- foreach(i=1:5, .options.RNG=123) %dorng% { runif(3) }
identical(s, s2)


## ----options_single_normalkind-------------------------------------------
## Pass the Normal RNG kind to use within the loop
# results are identical if not using the Normal kind in the loop
optsN <- list(123, normal.kind="Ahrens")
resN.U <- foreach(i=1:5, .options.RNG=optsN) %dorng% { runif(3) }
identical(resN.U[1:5], res[1:5])

# Results are different if the Normal kind is used and is not the same
resN <- foreach(i=1:5, .options.RNG=123) %dorng% { rnorm(3) }
resN1 <- foreach(i=1:5, .options.RNG=optsN) %dorng% { rnorm(3) }
resN2 <- foreach(i=1:5, .options.RNG=optsN) %dorng% { rnorm(3) } 
identical(resN[1:5], resN1[1:5])
identical(resN1[1:5], resN2[1:5])


## ----options_6length-----------------------------------------------------
# use a 6-length numeric
s <- foreach(i=1:5, .options.RNG=1:6) %dorng% { runif(3) }
attr(s, 'rng')[1:3]


## ----options_7length-----------------------------------------------------
# use a 7-length numeric, used as first value for .Random.seed
seed <- attr(res, 'rng')[[2]]
s <- foreach(i=1:5, .options.RNG=seed) %dorng% { runif(3) }
identical(s[1:4], res[2:5])


## ----options_list--------------------------------------------------------
# reproduce previous %dorng% loop
s <- foreach(i=1:5, .options.RNG=res) %dorng% { runif(3) }
identical(s, res)

## use completely custom sequence of seeds (e.g. using RNG "Marsaglia-Multicarry")
# as a matrix
seedM <- rbind(rep(401, 5), mapply(rep, 1:5, 2))
seedM
sM <- foreach(i=1:5, .options.RNG=seedM) %dorng% { runif(3) }
# same seeds passed as a list
seedL <- lapply(seq(ncol(seedM)), function(i) seedM[,i])
sL <- foreach(i=1:5, .options.RNG=seedL) %dorng% { runif(3) }
identical(sL, sM)


## ----set_seed_diff-------------------------------------------------------
# default RNG kind
RNGkind('default')
def <- foreach(i=1:5, .options.RNG=123) %dorng% { runif(3) }

# Marsaglia-Multicarry
RNGkind('Marsaglia')
mars <- foreach(i=1:5, .options.RNG=123) %dorng% { runif(3) }
identical(def, mars)

# revert to default RNG kind
RNGkind('default')


## ----schedule, tidy=FALSE------------------------------------------------
# define a stochastic task to perform
task <- function() c(pid=Sys.getpid(), val=runif(1))

# using the previously registered cluster with 2 workers
set.seed(123)
res_2workers <- foreach(i=1:5, .combine=rbind) %dorng% {
	task()
}
# stop cluster
stopCluster(cl)

# Sequential computation
registerDoSEQ()
set.seed(123)
res_seq <- foreach(i=1:5, .combine=rbind) %dorng% {
	task() 
}
#

# Using 3 workers
# NB: if re-running this vignette you should edit to force using 3 here 
cl <- makeCluster( if(isManualVignette()) 3 else 2)
length(cl)
# register new cluster
registerDoParallel(cl)
set.seed(123)
res_3workers <- foreach(i=1:5, .combine=rbind) %dorng% { 
	task()
}
# task schedule is different
pid <- rbind(res1=res_seq[,1], res_2workers[,1], res2=res_3workers[,1])
storage.mode(pid) <- 'integer'
pid
# results are identical
identical(res_seq[,2], res_2workers[,2]) && identical(res_2workers[,2], res_3workers[,2])


## ----registerDoRNG-------------------------------------------------------

set.seed(123)
res <- foreach(i=1:5) %dorng% { runif(3) }

registerDoRNG(123)
res_dopar <- foreach(i=1:5) %dopar% { runif(3) }
identical(res_dopar, res)
attr(res_dopar, 'rng')


## ----multiple, tidy=FALSE------------------------------------------------
set.seed(456)
s1 <- foreach(i=1:5) %dorng% { runif(3) }
s2 <- foreach(i=1:5) %dorng% { runif(3) }
# the two loops do not use the same streams: different results
identical(s1, s2)

# but the sequence of loops is reproducible as a whole
set.seed(456)
r1 <- foreach(i=1:5) %dorng% { runif(3) }
r2 <- foreach(i=1:5) %dorng% { runif(3) }
identical(r1, s1) && identical(r2, s2) 

# one can equivalently register the doRNG backend and use %dopar%
registerDoRNG(456)
r1 <- foreach(i=1:5) %dopar% { runif(3) }
r2 <- foreach(i=1:5) %dopar% { runif(3) }
identical(r1, s1) && identical(r2, s2)


## ----nested_error, error = TRUE, try = TRUE------------------------------
# nested loop
try( foreach(i=1:10) %:% foreach(j=1:i) %dorng% { rnorm(1) } )

# conditional loop
try( foreach(i=1:10) %:% when(i %% 2 == 0) %dorng% { rnorm(1) } )


## ----nested--------------------------------------------------------------

# doRNG must not be registered
registerDoParallel(cl)

# generate sequence of seeds of length the number of computations
n <- 10; p <- 5
rng <- RNGseq( n * p, 1234)

# run standard nested foreach loop
res <- foreach(i=1:n) %:% foreach(j=1:p, r=rng[(i-1)*p + 1:p]) %dopar% {
    
	# set RNG seed
    rngtools::setRNG(r)
	
    # do your own computation ...
    c(i, j, rnorm(1))
}

# Compare against the equivalent sequential computations
k <- 1
res2 <- foreach(i=1:n) %:% foreach(j=1:p) %do%{
    # set seed
	rngtools::setRNG(rng[[k]])
	k <- k + 1
	
    # do your own computation ...
	c(i, j, rnorm(1))
}

stopifnot( identical(res, res2) )


## ----nested_unequal------------------------------------------------------
# generate sequence of seeds of length the number of computations
n <- 10
rng <- RNGseq( n * (n+1) / 2, 1234)

# run standard nested foreach loop
res <- foreach(i=1:n) %:% foreach(j=1:i, r=rng[(i-1)*i/2 + 1:i]) %dopar%{
	
	# set RNG seed
	rngtools::setRNG(r)
	
	# do your own computation ...
	c(i, j, rnorm(1))
}

# Compare against the equivalent sequential computations
k <- 1
res2 <- foreach(i=1:n) %:% foreach(j=1:i) %do%{
	# set seed
	rngtools::setRNG(rng[[k]])
	k <- k + 1
	
	# do your own computation ...
	c(i, j, rnorm(1))
}

stopifnot( identical(res, res2) )


## ----conditional---------------------------------------------------------

# un-conditional single loop
resAll <- foreach(i=1:n, .options.RNG=1234) %dorng%{
	# do your own computation ...
	c(i, rnorm(1))
}

# generate sequence of RNG
rng <- RNGseq(n, 1234)

# conditional loop: even iterations
resEven <- foreach(i=1:n, r=rng) %:% when(i %% 2 == 0) %dopar%{
	
	# set RNG seed
	rngtools::setRNG(r)
	
	# do your own computation ...
	c(i, rnorm(1))
}

# conditional loop: odd iterations
resOdd <- foreach(i=1:n, r=rng) %:% when(i %% 2 == 1) %dopar%{
	
	# set RNG seed
	rngtools::setRNG(r)
	
	# do your own computation ...
	c(i, rnorm(1))
}

# conditional loop: only first 2 and last 2
resFL <- foreach(i=1:n, r=rng) %:% when(i %in% c(1,2,n-1,n)) %dopar%{
	
	# set RNG seed
	rngtools::setRNG(r)
	
	# do your own computation ...
	c(i, rnorm(1))
}

# compare results
stopifnot( identical(resAll[seq(2,n,by=2)], resEven) )
stopifnot( identical(resAll[seq(1,n,by=2)], resOdd) )
stopifnot( identical(resAll[c(1,2,n-1,n)], resFL) )



## ----nested_conditional--------------------------------------------------
# generate sequence of seeds of length the number of computations
n <- 10
rng <- RNGseq( n * (n+1) / 2, 1234)

# run standard nested foreach loop
res <- foreach(i=1:n) %:% when(i %% 2 == 0) %:% foreach(j=1:i, r=rng[(i-1)*i/2 + 1:i]) %dopar%{
	
	# set RNG seed
	rngtools::setRNG(r)
	
	# do your own computation ...
	c(i, j, rnorm(1))
}

# Compare against the equivalent sequential computations
k <- 1
resAll <- foreach(i=1:n) %:% foreach(j=1:i) %do%{
	# set seed
	rngtools::setRNG(rng[[k]])
	k <- k + 1
	
	# do your own computation ...
	c(i, j, rnorm(1))
}

stopifnot( identical(resAll[seq(2,n,by=2)], res) )


## ----perf, cache=TRUE----------------------------------------------------
# load rbenchmark
library(rbenchmark)

# comparison is done on sequential computations
registerDoSEQ()
rPar <- function(n, s=0){ foreach(i=1:n) %dopar% { Sys.sleep(s) } }
rRNG <- function(n, s=0){ foreach(i=1:n) %dorng% { Sys.sleep(s) } }

# run benchmark
cmp <- benchmark(rPar(10), rRNG(10)
			, rPar(25), rRNG(25)
			, rPar(50), rRNG(50)
			, rPar(50, .01), rRNG(50, .01)
            , rPar(10, .05), rRNG(10, .05)
			, replications=5)
# order by increasing elapsed time
cmp[order(cmp$elapsed), ]


## ----doRNGversion--------------------------------------------------------
doRNGversion('1.3')


## ----doRNGversion_revert-------------------------------------------------
doRNGversion(NULL)


## ----news, echo=FALSE, results='asis'------------------------------------
cat(paste(readLines(system.file('NEWS', package='doRNG')), collapse="\n"))


## ----stopCluster---------------------------------------------------------
stopCluster(cl)


## ----session_info, echo=FALSE, comment=NA--------------------------------
sessionInfo()


