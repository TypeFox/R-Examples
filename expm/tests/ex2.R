
#### Example matrices from the Matlab demos //  expAtv() examples

library(expm)

source(system.file("test-tools.R", package= "expm"), keep.source=FALSE)
## -> assertError()...

## --- 1 ---
## Here, all three {eigen; Taylor; Pade(Scaling & Squaring)} should do well

A1 <- rbind(0:2,
            c(0.5, 0, 1),
            2:0)
A1
ml1 <- lapply(c(4:10,20),
              function(order) expm(A1, "Pade", order=order))
for(k in seq_len(length(ml1) - 1))
    stopifnot(all.equal(ml1[[k]], ml1[[k + 1]], tolerance = 1e-12))

for(k in seq_len(length(ml1) - 1)) {
    print(all.equal(ml1[[k]], ml1[[k + 1]], tolerance = 0))
}

mA1 <- ml1[[4]]
stopifnot(all.equal(mA1,
                    matrix(c(5.3090812852106, 2.8087900904073, 5.1737460019740,
                             4.0012030182399, 2.8845155413486, 4.0012030182399,
                             5.5778402926177, 3.1930144369526, 5.7131755758543),
                           3, 3),
                           check.attributes = FALSE, tolerance = 1e-11))


## --- 2 ---
## Here, Taylor "fails":

## A matrix where the terms in the Taylor series become very large
## before they go to zero.
A2 <- rbind(c(-147, 72),
            c(-192, 93))
A2
(mA2 <- expm(A2, method="Pade"))
stopifnot(all.equal(mA2,
                    matrix(c(-0.099574136735459, -0.199148273470915,
                              0.074680602551593 , 0.149361205103183),
                           2, 2), check.attributes = FALSE, tolerance = 1e-11))
mA2.T <- expm(A2, method = "Taylor")
stopifnot(all.equal(mA2, mA2.T, tolerance=1e-10))
all.equal(mA2, mA2.T, tolerance=1e-14)#-> 3.2e-12  {MM: I think that's pretty good}

## --- 3 ---
## Here, Eigenvalues  must fail ("not a full set"):
A3 <- rbind(c(-1, 1),
            c(0, -1))
(mA3 <- expm(A3, method="Pade"))
assertError(expm(mA3, method="R_Eigen"))
em1 <- exp(-1)
stopifnot(all.equal(mA3, ## and the exact solution:
                    matrix(c(em1, 0, em1, em1), 2, 2),
                    check.attributes = FALSE, tolerance = 1e-14))

## using 'eps' instead of 0 :
## ---> see m2ex3() etc in ./exact-ex.R


## --- 4 ---
## Here, some version of do_expm() failed:
(m <- matrix(c(0,2:0),2))
## Eigenvalue decomposition:
d <- c(sqrt(2), -sqrt(2))
V <- rbind(c(sqrt(1/3), -sqrt(1/3)),
           c(sqrt(2/3),  sqrt(2/3)))
## ==>
IV <- rbind(c( sqrt(3/4), sqrt(3/8)),
            c(-sqrt(3/4), sqrt(3/8)))
stopifnot(all.equal(V %*% IV, diag(2)))
em.true <- V %*% (exp(d) * IV)
stopifnot(all.equal(em.true, expm::expm(m)),
          all.equal(em.true, expm::expm(m,"Pade"), check.attributes=FALSE))

###----------- expAtv() ----------------

## Bug report, 8 Sep 2014  (R-forge Bugs item #5919), by: Peter Ralph
stopifnot(expAtv(A3, v=c(0,0))$eAtv == 0)


n <- 500
A <- bandSparse(n,n, -1:1, diag = list(-(2:n), -5*(1:n), 1:(n-1)))
v <- 100*(n:1)
t.v <- showSys.time(rr <- expAtv(A, v=v))
if(doExtras) { ## this is an order of magnitude slower :
    t.A <- system.time(w. <- (eA <- expm(A, "Higham08")) %*% v)
    stopifnot(all.equal(rr$eAtv, as.numeric(w.)))
    print( mean((t.A / t.v)[c(1,3)]) )## 23.57 {nb-mm3}; 21.0 {lynne}
}


## Bug report on R-forge  by Peter Ralph (petrelharp)
## If the entries of A are less than about 1e-8, then expAtv(A,v) fails
## with Error: length(d <- dim(x)) == 2 is not TRUE
## ... an error that comes from expm, because it has got a 1x1 matrix. (I can't tell why this causes an error; outside of expAtv this doesn't cause an error...)

## To reproduce:

##' @title Compute several "scaled" versions of  e^{At} v :
##' @param A n x n matrix
##' @param v n  vector
##' @param s vector of scales
##' @return list of  expAtv() results
##' @author Martin Maechler, based on Peter Ralph's idea:
scl.e.Atv <- function(A, v, s) {
    c(list(I = expAtv(A, v)),
      unlist(lapply(s, function(l) {
          ## cat(sprintf(" %7g\n", l))
          list(lA = expAtv(l*A, v), lAI = expAtv(l*A, v, t=1/l))
      }), recursive = FALSE))
}

A <- matrix( 1:9, nrow=3 )/8
v <- rep(1,3)
sc <- 4^c(-500, -200, -100, -5*(15:6), -2*(14:9), -17:15)
## 10^9 is too large => expm() "overflow" NaN
r <- scl.e.Atv(A,v, s = sc) # at least without error
(eAv <- t(sapply(r, `[[`, "eAtv")))
## Ensure that indeed	expAtv(A, v)  =.=  expAtv(e*A, v, 1/e)  for e > 0
## -----                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eAv[1,]
assert.EQ.mat(unname( eAv[rep(1, length(sc)), ]),
	      unname( eAv[1+2*seq_along(sc), ] ), tolerance=1e-14)
					# 64-bit lynne: 2.7e-16 !!

sc.Atv <- function(A,v, s) {
    vapply(s, function(l) expAtv(l*A, v, t=1/l)$eAtv, v)
}

chk.sc.Atv <- function(A,v, s, tol=1e-15) {
    r <- vapply(s, function(l) expAtv(l*A, v, t=1/l)$eAtv, v)
    I <- expAtv(A,v)$eAtv
    if (!isTRUE(eq <- all.equal(as.vector(r), rep(I, length(s)), tolerance = tol)))
	stop("not all.equal() |->  ", eq)
}

chk.sc.Atv(A,v, sc, tol=1e-14)
## for information: see the precision:
tryCatch( chk.sc.Atv(A,v, sc, tol= 0), error=identity)$message


A0 <- matrix( c(-3,1,2,1,-2,1,0,1,-1), nrow=3, byrow=TRUE)
A1 <- A0 + 1e-16*rnorm(9)
## These two also failed originally
chk.sc.Atv(A0, v=10^(1:3), s=sc, tol=1e-14)
chk.sc.Atv(A1, v=rep(1,3), s=sc, tol=1e-14)

set.seed(17)
S <- rSpMatrix(29, density = 1/64)
v <- round(100*rnorm(nrow(S)))
if(FALSE) ## Error in  balance(baP$z, "S") :
    ## BLAS/LAPACK routine 'DGEBAL' gave error code -3
    chk.sc.Atv(S/64, v, s=sc, tol=1e-14)
if(FALSE) {
    ## after
    debug(chk.sc.Atv)
    ## this is revealing:
    image(as(relErrV(I, r),"sparseMatrix"))
    ## ==>
    sc[28:29] # are giving the largest differences
}
