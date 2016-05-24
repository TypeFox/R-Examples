require(simsalapar)

vl <- varlist(n.sim = list(type="N", expr = quote(N[sim]), value = 7),
              b     = list(type="grid", value = 1:3),
              a     = list(type="grid", value = 1:2),
              ii    = list(type="inner", value = 1:5))
(gr <- data.matrix(mkGrid(vl)))

do.one <- function(a,b, ii) {
    ## ii with names that are propagated (!) :
    names(ii) <- paste("I", ii, sep=".")
    t <- 10*(10*a + b) + round(runif(1), 1)/4
    ii + t
}

set.seed(17)
rL  <- doLapply(vl, doOne=do.one, repFirst=TRUE)
rL2 <- doLapply(vl, doOne=do.one, repFirst=FALSE)

va  <- getArray(rL , "value")
str(va)
va2 <- getArray(rL2, "value")
str(va2)

stopifnot( all.equal(va, va2, tol = .001) )


## approximate mean result
am <- outer(1:5, 10*outer(1:3, 10*(1:2), "+"), "+")

## Test 'repFirst=TRUE' :

m1 <- apply(va, 1:3, mean)
stopifnot(unname(dim(m1)) == dim(am)) # dim() match {apart from names}
dim     (am) <- dim     (m1)
dimnames(am) <- dimnames(m1)
stopifnot(all.equal(m1, am, tol = .001),
          apply(va, 1:3, sd) < 0.1)

## Test 'repFirst=FALSE' :

m2 <- apply(va2, 1:3, mean)
stopifnot(unname(dim(m2)) == dim(am)) # dim() match {apart from names}
stopifnot(all.equal(m2, am, tol = .001),
          apply(va2, 1:3, sd) < 0.1)

