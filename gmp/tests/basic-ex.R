library(gmp)

## From ~/R/Pkgs/Matrix/inst/test-tools-1.R:

##' @title Ensure evaluating 'expr' signals an error
##' @param expr
##' @return the caught error, invisibly
##' @author Martin Maechler
assertError <- function(expr) {
    d.expr <- deparse(substitute(expr))
    t.res <- tryCatch(expr, error = function(e) e)
    if(!inherits(t.res, "error"))
	stop(d.expr, "\n\t did not give an error", call. = FALSE)
    invisible(t.res)
}

Z1 <- as.bigz(1) ; Z1[FALSE]
Q1 <- as.bigq(1) ; Q1[FALSE]
stopifnot(0 == length(z0 <- as.bigz(0[FALSE])),# failed earlier
	  0 == length(q0 <- as.bigq(0[FALSE])),# ditto
	  is.bigz(Z1), is.bigz(z0), !is.bigz(1L), !is.bigz(1),   !is.bigz(Q1),
	  is.bigq(Q1), is.bigq(q0), !is.bigq(1L), !is.bigq(1/2), !is.bigq(Z1))

Z1[integer()] <- 2 # segfaulted earlier
Q1[integer()] <- 2 # ditto
assertError(Z1[1] <- list(1)) # segfaulted
assertError(Q1[1] <- list(1)) #  "
assertError(Z1[1] <- NULL   ) #  "
assertError(Q1[1] <- NULL   ) #  "

stopifnot(identical(Z1, as.bigz(1L)),  identical(Q1, as.bigq(1L)),
          identical(1L, as.integer(Z1)),
          identical(1L, as.integer(Q1)),## failed earlier
          identical(as.bigz(1[FALSE]), Z1[FALSE]),
          identical(as.bigz(1[-1]),    Z1[-1]),
          identical(Z1[-1], rep(Z1, 0))
          , ##----------- bigq -------------
          identical(as.bigq(1[FALSE]), Q1[-1]),
          identical(Q1[FALSE], Q1[-1]),
          identical(Q1[-1], rep(Q1, 0)),
          identical(q0, rep(Q1, 0))
          )

stopifnot(length(1[0]) == 0, 0 == length(Z1[0]))
Z <- as.bigz(I <- 2^(5*0:5)); mZ <- as.bigz(mI <- matrix(I, 2,3))
Q <- Z / 4                  ; mQ <- matrix(Q, 2,3)

ii <- c(3:2,0:2,1:0,0:2)
i. <- c(2:0,1:0,1); j. <- ii[1:7]
i <- i.[i. != 0]
j <- j.[j. != 0]
I[ii] ; mI[i.,j.]
stopifnot(all.equal(  Z[ii], I[ii], tol=0),
	  all.equal(4*Q[ii], I[ii], tol=0),
	  identical(mI[i,j], mI[i.,j.]),
	  identical(mZ[i,j], mZ[i.,j.]),
	  identical(mQ[i,j], mQ[i.,j.]))
stopifnot(all.equal(asNumeric(mZ[i,j]), mI[i,j], tol=0),
	  all.equal(	    4*mQ[i,j],	mI[i,j], tol=0))

## Outside indexing for *matrices* now gives an error:
assertError(mI[1,4]); assertError(mZ[1,4]); assertError(mQ[1,4])
assertError(mI[3,2]); assertError(mZ[3,2]); assertError(mQ[3,2])
## whereas outside indexing of vectors should give NA:
stopifnot(identical(I[8:5], asNumeric(Z[8:5])),
	  identical(I[8:5], asNumeric(Q[8:5] * 4)))

## "basics", including as.matrix(), as.array(), as.list() :
i <- 1:9
(x <- as.bigz(i, mod = 3))
mx <- as.matrix(x) ## used to "bomb" badly:
## (terminate called after throwing an instance of 'std::bad_alloc')
lx <- as.list(x)
stopifnot(5*x == (5*i) %% 3,
	  identical(mx, as.array(x)),
	  is(mx, "bigz"), dim(mx) == c(9,1),
	  is.list(lx),
	  identical(unlist(lx),
		    unlist(lapply(x, unclass))))

## remove modulus "the new way" (NULL did fail):
modulus(x) <- NULL
Q <- x / 2
mq <- as.matrix(Q)
lq <- as.list(Q)
stopifnot(identical(x, as.bigz(i %% 3)),
	  identical(mq, as.array(Q)),
	  is(mq, "bigq"), dim(mq) == c(9,1),
	  is.list(lq),
	  identical(unlist(lq),
		    unlist(lapply(Q, unclass))))

## duplicated(), unique() : ----------------------
q7 <- as.bigq(-5:7, 7)
if(FALSE)# not yet {well, *HARD* / impossible(?) without S4 }
Q <- q7^2 * as.bigz(77)^10
Q <- q7^2 * as.bigq(77, 2)^10
(uQ <- unique(Q))
(sDup <- sum(duplicated(Q))) # = 5
stopifnot(!duplicated(uQ),
          sDup + length(uQ) == length(Q))
nQ <- asNumeric(Q)

stopifnot( identical(duplicated(Q), duplicated(nQ))
	  , all.equal(unique(Q), unique(nQ))
          , sort(asNumeric(unique(denominator(Q)))) == 4^c(0, 3:5)
	  , TRUE)

## _ TODO _  rep()  [times, length.out, each]
checkRep <- function(x) {
    if((n <- length(x)) < 2) stop("'length(x)' must at least be 2, for these checks")
    ii <- seq_len(n)
    n1 <- pmin(.9*n, n-1)
    stopifnot(identical(rep(x, 1), x),
              identical(rep(x, 3), c(x,x,x)),
        identical(rep(x, length.out=n1), x[1:n1])
        ,
        identical(rep(x,    length.out=n+2), x[c(ii,1:2)])
        , ## times is *not* considered when 'length.out' is specified:
        identical(rep(x, 4, length.out=n+2), x[c(ii,1:2)])
        ,
        identical(rep(x, 2, length.out=n1), x[1:n1])
        ,
        identical(x, rep(x, each=2)[2*ii])
        )
}

checkRep(Q)
checkRep(q7)
(Nu <- numerator(uQ))
checkRep(Nu)

##------ Now check that  base :: pmin() / pmax()  works *in simple cases* for bigz
##------ (because  rep(., length.out) works:
## {{MM: compare with ~/R/Pkgs/Rmpfr/tests/arith-ex.R }}
(x <- as.bigz(ix <- 2^(3* 0:7)))
x9 <- pmin(x,9)
xp123 <- pmax(x, 123)
stopifnot(x9 == c(1,8, rep(9,6)),
          identical(x,  pmin(x, Inf)),
          identical(x9, pmin(x,23, Inf, 9)),
          xp123[1:3] == 123,
          xp123[-(1:3)] > 123)

if(FALSE)## FIXME
if(require("Rmpfr") && packageVersion("Rmpfr") >= "0.5-2") {
   stopifnot(
       all.equal(pmin(14,  x, 9),
                 pmin(14, ix, 9), tol=0)
       ,
       all.equal(mq <- pmin(14,  x/3, 9), ## numbers + bigq
                       pmin(14, ix/3, 9), tol= 1e-15)
       ,
       is.bigq(mq))
}
