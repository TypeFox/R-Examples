require("Rmpfr")
## includes ("gmp")# want to check "mixed arithmetic" too  __ TODO __

`%=N=%` <- function(x,y) (x == y) | (is.na(x) & is.na(y))
all.EQ <- function(x,y, tolerance = 2^-98, ...) # very small tol. for MPFR
    all.equal(x, y, tolerance=tolerance, ...)
warningI <- function(...) warning(..., immediate. = TRUE)

## Check that we got the "which.*" methods also from "bigq":
bcl <- c("ANY", "bigq", "bigz", "mpfr")
##if(packageVersion("gmp") >= "0.5-8") {
stopifnot(identical(bcl,
		    sort(unlist(findMethods("which.max")@signatures))),
	  identical(bcl,
		    sort(unlist(findMethods("which.min")@signatures))))
##}

options(warn = 1)# warnings *immediately*
(doExtras <- Rmpfr:::doExtras())
eps2 <-   2 * .Machine$double.eps
eps8 <-   8 * .Machine$double.eps
eps32 <- 32 * .Machine$double.eps

## must take the *larger* of the two precisions:
stopifnot(format(mpfr(1, 60) / mpfr(7, 160)) ==
          "0.14285714285714285714285714285714285714285714285712")

(x <- mpfr(0:7, 100) / 7)
ix <- x^-1000
iX <- asNumeric(ix)

stopifnot( mpfrIs0(x - x), # badly failed on 64-bit
	  identical(-x, 0-x),# testing "- x"
          all.equal(ix, (1/x)^1000, tol= 1e-25),
          is.numeric(iX), iX[1:4] == Inf, # failed previously as we used RNDD (downward rounding)
          all.equal(log(iX[5:8]), c(559.6157879, 336.4722366, 154.1506798, 0),
                    tol = 1e-9))

## checking hexadecimal input :
stopifnot(mpfr("0xFFFFFFFFFFFFFFFFFFFF", base=16) + 1 == 2^80,
## sign(0) == 0:
          identical(sign(as(-1:1, "mpfr")), -1:1 + 0))

stopifnot(all.equal(as.numeric(x+ 1L),
                    as.numeric(x)+1L, tol = eps2),
	  as.integer(  x [x < 1]) == 0,# was *wrong* {we round()ed; previously "down"!}
	  as.integer((-x)[x < 1]) == 0,#  (ditto)
          (3 * x)/3 <= x,
          all.equal(as.numeric(x * 2L),
                    as.numeric(x + x), tol = 0))

u <- mpfr(0:17, 128)/17
two <- mpfr(2,100)
stopifnot(all.EQ(u ^ two, u ^ 2),
          identical(u ^ 2, u ^ 2L),
          all.EQ(two ^ u, 2 ^ u),
          identical(2 ^ u, 2L ^ u),
          floor  (3*u) == floor  (3/17*(0:17)),
          ceiling(u*5) == ceiling(5/17*(0:17))
          )

i7  <- mpfr(0:7,  200)/ 7
i17 <- mpfr(0:17, 300)/17
stopifnot(all.equal(as.numeric(x+1),
		    as.numeric(x)+1),
	  all.equal(round(x,2), round(asNumeric(x), 2), tol=1e-15),
	  all.equal(round(mpfr(1.152, 80), 2), 1.15), # was wrong {as.integer() bug}
	  all.equal(0:7,   7 * round ( i7, 25), tol = 2e-25),
	  all.equal(0:7,   7 * round ( i7, 50), tol = 2e-50),
	  all.equal(0:17, 17 * signif(i17,100), tol = 2e-100),
	  all.equal(0:17, 17 * signif(i17, 20), tol = 2e-20)
          )

## When we compute with 100 bits,
## we should compare relative errors with  2^-100 :
del <- abs((x+pi)-pi - x) / 2^-100
stopifnot(del <= 4) ## <= 2 already
(fd <- format(del, drop0 = TRUE))
if(print(Sys.info()[["machine"]]) == "x86_64")
    stopifnot(fd %in% as.character(c(0:2, c(2,7)/4)))


checkPmin <- function(x, nx = as(x, "numeric")) {
    rx <- if(is(x,"mpfr")) round(x, 25) else x
    isZ <- is(x, "bigz") || is(nx, "bigz")
    M.X <- max(x, na.rm=TRUE)
    m.x <- min(x, na.rm=TRUE)
    stopifnot(all.equal(x, nx),
	      pmin(x, x, M.X) %=N=% x, x %=N=% pmax(x, m.x, x),
	      all.equal(x, pmin(x, nx, x, M.X)),
	      all.equal(x, pmax(m.x, nx, x, rx, m.x)),
	      if(isZ)TRUE else all.equal(pmin(x, 0.75), pmin(nx, 0.75)),
	      if(isZ)TRUE else all.equal(pmax(x, 0.25), pmax(nx, 0.25)))
}

x <- mpfr(0:7, 100) / 7
checkPmin(x)

nx <- (0:7)/7
(qx <- as.bigq(0:7, 7))
 x[c(2,5)] <- NA
nx[c(2,5)] <- NA
qx[c(2,5)] <- NA

Z <- as.bigz(1:7)
mZ <- mpfr(Z, 64)
stopifnot(Z == mZ, mZ == Z)

checkPmin(x, nx)
    cat("checking pmin(. bigq ): ")
    ## FIXME checkPmin(x, qx); cat("[Ok]\n")
    ##
    print( base::pmin(Z, Z, max(Z)) )# via  gmp:: rep.bigz(x, length.out = *)
    cat("checking pmin(. bigz ) : ")
    checkPmin(Z); cat("[Ok]\n") # via gmp:: all.equal.bigz()

stopifnot(all.equal( round(x, 10),  round(nx, 10)),
          all.equal(signif(x, 10), signif(nx, 10)))

## L & x ,  x & L  failed in Rmpfr 0.2* and 0.4-2
stopifnot(identical(L <- x > 0.5, L & x),
	  identical(L, x & L),
	  identical(x > 0, x | L))

##-------------- Modulo and "integer division" -------------

## R's	?Arithmetic :
##
##   ‘%%’ indicates ‘x mod y’ and ‘%/%’ indicates integer division.  It
##   is guaranteed that ‘x == (x %% y) + y * ( x %/% y )’ (up to
##   rounding error) unless ‘y == 0’ where the result of ‘%%’ is
##   ‘NA_integer_’ or ‘NaN’ (depending on the ‘typeof’ of the
##   arguments).
##
## and has 'details' about  how	 non-integer 'y' works
##
N <- if(doExtras) 1000 else 200
mm <- c(-4:4, sample(50, N-9, replace=TRUE))
for(n in seq_len(N)) {
    cat("."); if(n %% 50 == 0) cat(n,"\n")
    m <- mm[n]
    prec <- sample(52:200, 1)# "high precision" ==> can use small tol
    x <- sample(100, 50) - 20
    for(kind in c('int','real')) {
	if(kind == "real") {
	    m <- jitter(m)
	    x <- jitter(x)
	    tol.1 <- eps32 * pmax(1, 1/abs(m))
	    EQ <- function(x,y, tol = tol.1)
		isTRUE(all.equal(x, as.numeric(y), tol=tol))
	    EQ2 <- function(x,y, tol = tol.1) {
		## for the DIV - MOD identity, a small x leads to cancellation
		all((x %=N=% y) | abs(x - y) < tol*pmax(abs(x), 1)) ||
		    isTRUE(all.equal(x, as.numeric(y), tol=tol))
	    }
	} else { ## "integer"
	    EQ2 <- EQ <- function(x,y, tol) all(x %=N=% y)
	}
	i.m <- mpfr(x, prec) %% mpfr(m, prec)
	if(!EQ2(x %% m, i.m)) {
	    cat("\n -- m = ",m,"  (prec = ",prec,")\n")
	    rE <- range(rel.E <- as.numeric(1 - (x %% m)/i.m))
	    print(cbind(x, 'R.%%' = x %% m, rel.E))
	    MSG <- if(max(abs(rE)) < 1e-10) warningI else stop
	    MSG(sprintf("not all equal: range(rel.Err.) = [%g, %g]", rE[1],rE[2]))
	}
	##
	if(m != 0) {
	    ##---Check the    x	 ==  (x %% m) +	 m * ( x %/% m )  assertion ------
	    ##
	    if(EQ2(x, (x %% m) + m*( x %/% m ), tol = 1e-12)) { ## ok for R
		## --> also ok for mpfr ?
		iDm <- mpfr(x, prec) %/% mpfr(m, prec)
		rhs <- i.m + m*iDm
		if(!EQ2(x, i.m + m*iDm)) {
		    cat("\n -- m = ",m,"  (prec = ",prec,")\n")
		    print(cbind(x,' MPFR[ x%%m + m(x %/% m) ]' = as.numeric(rhs), rel.E))
		    MSG <- if(max(abs(rE)) < 1e-10) warningI else stop
		    MSG(sprintf("Identity(MOD - DIV) not all eq.: range(rel.Err.) = [%g, %g]",
				rE[1],rE[2]))
		}
	    } else {
		cat("\n hmm.. the basic	 %% <-> %/%  assertion 'fails' in *R* :\n")
		rhs <- (x %% m)	 +  m * ( x %/% m )
		rel.E <- (1 - rhs/x)
		print(cbind(x, 'x%%m + m(x %/% m)' = rhs, rel.E))
	    }
	}
    }
}

## mpfr  o  <number>  now implemented, for  '%%', too :
r <- as.double(i <- -10:20)

stopifnot(
    ## %% -------------------------------------
    mpfr(i, prec=99) %% 7  == i %% 7
    , ##
    mpfr(i, prec=99) %% 7  ==
    mpfr(i, prec=99) %% 7L
    , ##
    i %% mpfr(27, prec=99) == i %% 27
    , ##
    r %% mpfr(27, prec=99) == r %% 27
    , ## %/% -------------------------------------
    mpfr(i, prec=99) %/% 7  == i %/% 7
    , ##
    mpfr(i, prec=99) %/% 7  ==
    mpfr(i, prec=99) %/% 7L
    , ##
    mpfr(i, prec=99) %/% mpfr(27, prec=99) == i %/% 27
    , ##
    i %/% mpfr(27, prec=99) == i %/% 27
    , ##
    i %/% mpfr(27, prec=99) ==
    r %/% mpfr(27, prec=99)
    , TRUE ##
)

cat('Time elapsed: ', proc.time(),'\n') # "stats"

###------Standard Statistics Functions --------------------------------------------------------

x <- c(del, 1000)
stopifnot(identical(mean(x), mean(x, trim=0)))
for(tr in (0:8)/16)
    stopifnot(all.equal(mean(          x,  trim = tr),
                        mean(asNumeric(x), trim = tr), tol=1e-15))

cat('Time elapsed: ', proc.time(),'\n') # "stats"
