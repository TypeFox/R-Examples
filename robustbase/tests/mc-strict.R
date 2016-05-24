
#### Testing  medcouple	 and related functions

### here, we do "strict tests" -- hence no *.Rout.save
### hence, can also produce non-reproducible output such as timing

library(robustbase)
source(system.file("xtraR/mcnaive.R", package = "robustbase"))# mcNaive()
source(system.file("test-tools-1.R",  package="Matrix", mustWork=TRUE))
assertEQm12 <- function(x,y, giveRE=TRUE, ...)
    assert.EQ(x,y, tol = 1e-12, giveRE=giveRE, ...)
## ^^ shows *any* difference ("tol = 0") unless there is no difference at all
##
c.time <- function(...) cat('Time elapsed: ', ..., '\n')
S.time <- function(expr) c.time(system.time(expr))
DO <- function(...) S.time(stopifnot(...))

n.set <- c(1:99, 1e5L+ 0:1) # large n gave integer overflow in earlier versions
DO(0 == sapply(n.set, function(n) mc(seq_len(n))))
DO(0 == sapply(n.set, function(n) mc(seq_len(n), doRefl=FALSE)))

DO(0 == sapply(1:100, function(n) mcNaive(seq_len(n), "simple")))
DO(0 == sapply(1:100, function(n) mcNaive(seq_len(n), "h.use" )))


x1 <- c(1, 2, 7, 9, 10)
mcNaive(x1) # = -1/3
assertEQm12(-1/3, mcNaive(x1))
assertEQm12(-1/3, mcNaive(x1, "h.use"))
assertEQm12(-1/3, mc(x1))

x2 <- c(-1, 0, 0, 0, 1, 2)
mcNaive(x2, meth="simple") # = 0 - which is wrong
mcNaive(x2, meth="h.use")  # = 1/6 = 0.16666
assertEQm12(1/6, mc(x2))
assertEQm12(1/6, mcNaive(x2, "h.use"))

x4 <- c(1:5,7,10,15,25, 1e15) ## - bombed in orignal algo
mcNaive(x4,"h.use") # 0.5833333
assertEQm12( 7/12, mcNaive(x4, "h.use"))
assertEQm12( 7/12, mc( x4, doRefl= FALSE))
assertEQm12(-7/12, mc(-x4, doRefl= FALSE))


set.seed(17)
for(n in 3:50) {
    cat(" ")
    for(k in 1:5) {
	x <- rlnorm(n)
	mc1 <- mc(x)
	mc2 <- mcNaive(x, method = "simple")
	mc3 <- mcNaive(x, method = "h.use" )
	stopifnot(all.equal(mc1, mc3, tolerance = 1e-10),# 1e-12 not quite ok
		  mc2 == mc3)
	cat(".")
    }
};  cat("\n")

###----  Strict tests of adjOutlyingness():
###                      ================= changed after long-standing bug fix in Oct.2014

set.seed(1);  S.time(a1.1 <- adjOutlyingness(longley))
set.seed(11); S.time(a1.2 <- adjOutlyingness(longley))
##
set.seed(2); S.time(a2 <- adjOutlyingness(hbk))
set.seed(3); S.time(a3 <- adjOutlyingness(hbk[, 1:3]))# the 'X' space
set.seed(4); S.time(a4 <- adjOutlyingness(milk)) # obs.63 = obs.64
set.seed(5); S.time(a5 <- adjOutlyingness(wood))
set.seed(6); S.time(a6 <- adjOutlyingness(wood[, 1:5]))# the 'X' space

## 32-bit <-> 64-bit different results {tested on Linux only}
is32 <- .Machine$sizeof.pointer == 4 ## <- should work for Linux/MacOS/Windows
isMac <- Sys.info()["sysname"] == "Darwin"
isSun <- Sys.info()["sysname"] == "SunOS"
Rnk <- function(u) rank(unname(u), ties.method = "first")
## to use for testing below:
cat("\nRnk(a3 $ adjout): "); dput(Rnk(a3$adjout), control= {})
cat("\nRnk(a4 $ adjout): "); dput(Rnk(a4$adjout), control= {})

stopifnot(which(!a2$nonOut) == 1:14,
	  which(!a3$nonOut) == 1:14,
	  if(isSun || isMac || is32) TRUE else
	  ## which(!a4$nonOut) == if(is32 && !isMac) c(1, 2, 41, 70) else c(12, 70),
          which(!a4$nonOut) == c(9:19, 23:27,57, 59, 70, 77),
	  ## 'longley', 'wood' have no outliers in the "adjOut" sense:
	  ## FIXME: longley is platform dependent too
	  if(isMac) TRUE else sum(a1.2$nonOut) >= 15, # sum(.) = 16 [nb-mm3, Oct.2014]
	  a5$nonOut,
          a6$nonOut[-20],
	  ## hbk (n = 75) :
	  abs(Rnk(a3$adjout) -
             c(62, 64, 68, 71, 70,   65, 66, 63, 69, 67,   73, 75, 72, 74, 25,
               52, 44,  5, 11, 33,    6, 21, 29, 28, 59,    9, 12, 13, 37, 27,
               43, 35, 22, 55, 14,    2, 26, 46, 54, 15,   23, 41, 40, 32, 60,
               30, 61, 19, 16,  8,   39, 53, 51, 48, 20,   47, 50, 42,  7, 38,
               17, 57, 45, 18, 24,   34,  3, 58, 56,  4,    1, 10, 31, 36, 49)
	      ) <= 3 ## all 0 on 32-bit Linux
         ,
	  ## milk (n = 86) : -- Quite platform dependent!
      {
	  r <- Rnk(a4$adjout)
	  r64 <- ## the 64-bit (ubuntu 14.04, nb-mm3) values:
	      c(65, 66, 61, 56, 47,   51, 19, 37, 74, 67,   79, 86, 83, 84, 85,
		82, 81, 73, 80, 55,   27,  3, 70, 68, 78,   76, 77, 53, 48,  8,
		29, 33,	 6, 32, 28,   31, 36, 40, 22, 58,   64, 52, 39, 63, 44,
		30, 57, 46, 43, 45,   25, 54, 12,  1,  9,    2, 71, 14, 75, 23,
		 4, 10, 34, 35, 17,   24, 15, 20, 38, 72,   42, 13, 50, 60, 62,
		26, 69, 18,  5, 21,    7, 49, 11, 41, 59,   16)
          r32 <- ## Linux 32bit (florence: 3.14.8-100.fc19.i686.PAE)
              c(78, 79, 72, 66, 52,   61, 22, 41, 53, 14,   74, 85, 82, 83, 84,
                80, 81, 56, 73, 65,   30,  3, 16, 17, 68,   57, 58, 63, 54,  8,
                32, 37,  6, 36, 31,   35, 40, 44, 25, 69,   77, 62, 43, 76, 48,
                34, 67, 51, 47, 49,   28, 64, 12,  1,  9,    2, 33, 15, 59, 26,
                 4, 10, 38, 39, 20,   27, 18, 23, 42, 86,   46, 13, 60, 71, 75,
                29, 50, 21,  5, 24,    7, 55, 11, 45, 70,   19)
          d <- (r - if (is32) r32 else r64)
          if(has.d <- any(d != 0)) { print(cbind(r, d)); print(table(abs(d))) }
	  ## for the biggest part (79 out of 86), the ranks are "close":
          ## 2014: still true, but in a different sense..
          if(has.d)
              sum(abs(d) <= 17) >= 78 && sum(abs(d) <= 13) >= 75
          else TRUE
      })


## check of adjOutlyingness *free* bug
## reported by Kaveh Vakili <Kaveh.Vakili@wis.kuleuven.be>
set.seed(-37665251)
X <- matrix(rnorm(100*5),100,5)
Z <- matrix(rnorm(100*5,0,1/100),10,5)
Z <- sweep(Z, 2, c(5,rep(0,4)), FUN="+")
X[91:100,] <- Z
for (i in 1:10) {
    ## this would produce an error in the 6th iteration
    aa <- adjOutlyingness(x=X,ndir=250)
}

## "large n" (this did overflow sum_p, sum_q  earlier ==> had inf.loop):
set.seed(3); x <- rnorm(2e5)
(mx <- mc(x, trace.lev=3))
stopifnot(print(abs(mx - -0.000772315846101988)) < 1e-15)
					# 3.252e-19, 64b Linux
					# 1.198e-16, 32b Windows

### Some platform info :
local({ nms <- names(Si <- Sys.info())
        dropNms <- c("nodename", "machine", "login")
        structure(Si[c("nodename", nms[is.na(match(nms, dropNms))])],
                  class="simple.list") })

if(identical(1L, grep("linux", R.version[["os"]]))) { ##----- Linux - only ----
    ##
    Sys.procinfo <- function(procfile)
    {
        l2 <- strsplit(readLines(procfile),"[ \t]*:[ \t]*")
        r <- sapply(l2[sapply(l2, length) == 2],
                    function(c2)structure(c2[2], names= c2[1]))
        attr(r,"Name") <- procfile
        class(r) <- "simple.list"
        r
    }
    ##
    Scpu <- Sys.procinfo("/proc/cpuinfo")
    Smem <- Sys.procinfo("/proc/meminfo")
    print(Scpu[c("model name", "cpu MHz", "cache size", "bogomips")])
    print(Smem[c("MemTotal", "SwapTotal")])
}

##' Checking the breakdown point of mc() --- Hubert et al. theory said : 25%
##' using non-default  doReflect=FALSE  as that corresponds to original Hubert et al.
##'
##' @title Medcouple mc() checking
##' @param x
##' @param Xfun
##' @param eps
##' @param NAiferror
##' @param doReflect
##' @param ...
##' @return mc(*,..) or NaN in case mc() signals an error [non-convergence]
##' @author Martin Maechler
mcX <- function(x, Xfun, eps=0, NAiferror=FALSE, doReflect=FALSE, ...) {
    stopifnot(is.numeric(x), is.function(Xfun), "eps" %in% names(formals(Xfun)))
    myFun <-
	if(NAiferror)
	    function(u) tryCatch(mc(Xfun(u, eps=eps), doReflect=doReflect, ...),
				 error = function(e) NaN)
	else
	    function(u) mc(Xfun(u, eps=eps), doReflect=doReflect, ...)
    vapply(x, myFun, 1.)
}

X1. <- function(u, eps=0) c(1,2,3, 7+(-10:10)*eps, u + (-1:1)*eps)
## ==> This *does* breakdown [but points are not "in general position"]:
r.mc1 <- curve(mcX(x, X1.), 10, 1e35, log="x", n=1001)
rt1 <- uniroot(function(x) mcX(exp(x), X1.) - 1/2, lower=0, upper=500)
exp(rt1$root) #  4.056265e+31

## eps > 0  ==> No duplicated points ==> theory says breakdown point = 0.25
## -------  but get big numerical problems:
if(FALSE) { # ==> convergence problem [also in maxit = 1e5] .. really an *inf* loop!
r.mc1.1  <- curve(mcX(x, X1., eps= .1  ), 10, 1e35, log="x", n=1001)
r.mc1.2  <- curve(mcX(x, X1., eps= .01 ), 10, 1e35, log="x", n=1001)
r.mc1.3  <- curve(mcX(x, X1., eps= .001), 10, 1e35, log="x", n=1001)
r.mc1.5  <- curve(mcX(x, X1., eps= 1e-5), 10, 1e35, log="x", n=1001)
r.mc1.8  <- curve(mcX(x, X1., eps= 1e-8), 10, 1e35, log="x", n=1001)
r.mc1.15 <- curve(mcX(x, X1., eps=1e-15), 10, 1e35, log="x", n=1001)# still!
}
## practically identical to  eps = 0 where we have breakdown (see above)
r.mc1.16 <- curve(mcX(x, X1., eps=1e-16), 10, 1e35, log="x", n=1001)
all.equal(r.mc1, r.mc1.16, tol=1e-15)#-> TRUE

## Quite bad case: Non convergence
X2. <- function(u) c(1:3, seq(6, 8, by = 1/8), u, u, u)
try(mc(X2.(4.3e31)))## -> error: no convergence
if(FALSE) # and the same here -- after longer waiting:
    mc(X2.(4.3e31), eps1=1e-7, eps2=1e-100, maxit = 1e6)## -> error: no convergence

## related, more direct:
X3. <- function(u) c(10*(1:3), 60:80, (4:6)*u)
mc(X3.(1e31), trace=5) # fine convergence in one iter.
try(
mc(X3.(1e32), trace=3) # no convergence...
)# bad

try(mc(X3.(1e32), trace=5, maxit=6)) # no convergence...

### TODO : find example with *smaller* sample size -- with no convergence
X4. <- function(u, eps, ...) c(10, 70:75, (2:3)*u)
mc(X4.(1e34))# "fine"
## whoa: jump down and up:
r.mc4 <- curve(mcX(x, X4.), 100, 1e35, log="x", n=2^12)

X5. <- function(u) c(10*(1:3), 70:78, (4:6)*u)
try(mc(X5.(1e32), maxit=1000))

X5. <- function(u, eps,...) c(5*(1:12), (4:6)*u)
(r.mc5 <- mc(X5.(1e32), doReflect=FALSE, maxit=1000))
all.equal(1, ## <- i.e. complete breakdown
          r.mc5) ## platform dependent! yes, on 64-bit
try(mc(X5.(5e31), maxit=10000)) # no convergence..
r.mc5Sml <- curve(mcX(x, X5.), 1,  100, log="x", n=1024) ## quite astonishing
r.mc5Lrg <- curve(mcX(x, X5.), 1, 1e30, log="x", n=1024) ## ok..
## but then going higher -- we have problems:
r.mc5Big <- curve(mcX(x, X5., NAiferror=TRUE), 1, 1e38, log="x",
                  n = 2^12, type = "o", cex = 1/4)
warnings()
summary(r.mc5Big$y)
## 15 NA's at x :
with(r.mc5Big, x[is.na(y)])
## ~= [4.3, 5.8] * 10^31


c.time(proc.time())
