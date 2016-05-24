library(robustbase)

source(system.file("xtraR/test_MCD.R", package = "robustbase"))#-> doMCDdata
##          ../inst/xtraR/test_MCD.R
source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
## -> assertError(), relErr(), and:
showProc.time()

## -- now do it:
options(digits = 5)
set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed
doMCDdata()
doMCDdata(method="DetMCD"); warnings()
##                        vvvv no timing for 'R CMD Rdiff' outputs
doMCDdata(nrep = 12, time=FALSE)
doMCDdata(nrep = 12, time=FALSE, method="DetMCD"); warnings()
doMCDdata(nrep = 12, time=FALSE, method = "MASS")

###--- now the "close to singular" mahalanobis case:
set.seed(6)
(c3  <- covMcd(mort3))
(c3. <- covMcd(mort3, nsamp="deterministic"))
stopifnot(log(c3$crit) <= log(c3.$crit),
          print(log(c3.$crit / c3$crit)) <= 0.8)
## see 0.516 / 0.291 {with seed 7}
##
## rescale variables:
scaleV <- c(0.1, 0.1, 1, 1, .001, 0.1, 0.1, 100)
mm <- data.matrix(mort3) * rep(scaleV, each = nrow(mort3))
C3  <- covMcd(mm)
C3. <- covMcd(mm, nsamp="deterministic")
stopifnot(C3$mcd.wt == c3$mcd.wt)# here, not for all seeds!

## error ("computationally singular") with old (too high) default tolerance:
try( covMcd(mm, control= rrcov.control(tol = 1e-10)) )
try( covMcd(mm, control= rrcov.control(tol = 1e-10), nsamp="deterministic") )

showProc.time()

## "large" examples using different algo branches {seg.fault in version 0.4-4}:

n <- 600 ## - partitioning will be triggered
set.seed(1)
X <- matrix(round(100*rnorm(n * 3)), n, 3)
(cX  <- covMcd(X))
 cX. <- covMcd(X, nsamp="deterministic", scalefn = scaleTau2)
i <- names(cX); i <- i[!(i %in% c("call", "nsamp", "method", "raw.weights"))]
stopifnot(sum(cX.$raw.weights != cX$raw.weights) <= 2,
          all.equal(cX[i], cX.[i], tol= 1/9))

n <- 2000 ## - nesting will be triggered
set.seed(4)
X <- matrix(round(100*rnorm(n * 3)), n, 3)
set.seed(1)
summary(cX  <- covMcd(X)) # <- show newly activated  print.summary.mcd(.)
 cX. <- covMcd(X, nsamp="deterministic", scalefn = scaleTau2)
i2 <- i[i != "mcd.wt"]
stopifnot(print(sum(cX.$raw.weights != cX$raw.weights)) <= 3, # 2
          all.equal(cX[i2], cX.[i2], tol= 1/10))# 1/16

set.seed(1) ## testing of 'raw.only' :
cXo <- covMcd(X, raw.only=TRUE)
i <- paste0("raw.", c("cov", "center", "cnp2"))
stopifnot(cXo$raw.only, all.equal(cX[i], cXo[i], tol = 1e-15),
          c("best", "mah") %in% setdiff(names(cX), names(cXo)))
showProc.time()

## Now, some small sample cases:

## maximal values:
n. <- 10
p. <-  8
set.seed(44)
(X. <- cbind(1:n., round(10*rt(n.,3)), round(10*rt(n.,2)),
             matrix(round(10*rnorm(n. * (p.-3)), 1),  nrow = n., ncol = p.-3)))

## 2 x 1 ---> Error
r <- tryCatch(covMcd(X.[1:2, 2, drop=FALSE]), error=function(e)e)
stopifnot(inherits(r, "error"),
          grepl("too small sample size", r$message))

## 3 x 2 --- ditto
r <- tryCatch(covMcd(X.[1:3, 2:3]), error=function(e)e)
stopifnot(inherits(r, "error"),
          grepl("too small sample size", r$message))

## 5 x 3  [ n < 2 p  ! ]  --- also works for MASS
X <- X.[1:5, 1:3]
set.seed(101)
## the finite-sample correction is definitely doubtful:
summary(cc <- covMcd(X, use.correction = FALSE))
str(cc) ## best = 2 3 4 5
mcc <- MASS::cov.mcd(X)
stopifnot(cc$best == mcc$best,
          all.equal(cc$center, mcc$center, tolerance = 1e-10),
          all.equal(c(mcc$cov / cc$raw.cov), rep(0.673549282206, 3*3)))

## p = 4 -- 6 x 4 & 7 x 4  [ n < 2 p  ! ]
p <- 4
n <- 7
X <- X.[1:n, 1+(1:p)]
stopifnot(dim(X) == c(n,p))
(cc <- covMcd(X, use.correction = FALSE))
str(cc) ## best = 1 2 4 5 6 7
mcc <- MASS::cov.mcd(X)
stopifnot(cc$best == mcc$best,
          all.equal(cc$center, mcc$center, tolerance = 1e-10),
          all.equal(c(mcc$cov / cc$raw.cov), rep(0.7782486992881, p*p)))
n <- 6
X <- X[1:n,]
(cc <- covMcd(X, use.correction = FALSE))
mcc <- MASS::cov.mcd(X)
stopifnot(cc$best == mcc$best,
          all.equal(cc$center, mcc$center, tolerance = 1e-10),
          all.equal(c(mcc$cov / cc$raw.cov), rep(0.7528695976179, p*p)))

showProc.time()

## nsamp = "exact" -- here for p=7
coleman.x <- data.matrix(coleman[, 1:6])
showSys.time(CcX <- covMcd(coleman.x, nsamp= "exact"))
showSys.time(Ccd <- covMcd(coleman.x, nsamp= "deterministic"))
stopifnot(all.equal(CcX$best,
		    c(2, 5:9, 11,13, 14:16, 19:20), tolerance=0),
	  intersect(CcX$best, Ccd$best) == c(2,5,7,8,13,14,16,19,20),
          relErr(CcX$crit, Ccd$crit) < 0.35 # see ~ 0.34
)
summary(Ccd)


demo(determinMCD)## ../demo/determinMCD.R
##   ----------- including simple "exactfit" (code = 3)
warnings()

showProc.time()
if(!robustbase:::doExtras()) quit()

## if ( doExtras ) -----------------------------------------------------------------
## ==============

##  (nmini, kmini) examples:
set.seed(7) ; X1 <- gendata(10000, p=13, eps = 0.30)
showSys.time(c1 <- covMcd(X1$X)) # 0.87 sec
chk.covMcd <- function(ans, ind) {
    stopifnot(inherits(ans, "mcd"))
    ## check that all outliers were detected:
    mod.outl <- which(ans$mcd.wt == 0)
    outl.detected <- (ind %in% mod.outl)
    if(!all(outl.detected)) {
        cat("The following outliers are *not* detected:\n")
        print(which(!outl.detected))
    }
    fp <- !(mod.outl %in% ind)
    if(any(fp)) {
        cat(sprintf("False positive \"outliers\" (a few expected) %d of %d (= %.2f%%):\n",
       	     sum(fp), nobs(ans), 100*sum(fp)/nobs(ans)))
        print(which(fp))
    } else cat("** No ** false positive outliers -- exceptional!\n")
}
##
chk.covMcd(c1, X1$xind)
cat("\ncovMcd(*, kmini=12, trace=2) ...\n------\n")
showSys.time(c2 <- covMcd(X1$X, kmini=12, trace=2))# slower..
chk.covMcd(c2, X1$xind)
## Comparing:
ii <- !(names(c1) %in% c("call", "method"))
cat("\ncovMcd(*, nsamp=\"deterministic\")\n")
showSys.time(cD <- covMcd(X1$X, nsamp="deterministic"))# quite slower than FASTMCD
chk.covMcd(cD, X1$xind)
cat("<.>$crit = log(det(.)) [minimal = best] :\n")
print(cbind(sort(c(default = c1$crit, kmini.12 = c2$crit, determin = cD$crit))))
i2 <- names(c1)[ii]; i2 <- i2[i2 != "nsamp"]
## closer coincidence if "raw.*" are dropped:
i3 <- i2; i3 <- i3[ - grep("^raw", i3) ]
stopifnot(all.equal(c1[ii], c2[ii], tol= 0.02),
          all.equal(cD[i2], c1[i2], tol= 0.02),
          all.equal(cD[i3], c1[i3], tol= 6e-4), # 4.60e-4
          ## the 0/1 weights coincide :
          cD$mcd.wt == c1$mcd.wt,
          c2$mcd.wt == c1$mcd.wt)
showProc.time()

## Radarexample --- already some in  ../man/radarImage.Rd <<<-------------
data(radarImage)
print(d <- dim(radarImage)); n.rI <- d[1]
## The 8 "clear" outliers (see also below)
ii8 <- c(1548:1549, 1553:1554, 1565:1566, 1570:1571)
set.seed(7)
showSys.time( L1 <- lapply(0:200, function(n)
    n+ which(0 == covMcd(unname(radarImage[(n+1L):n.rI,]), trace=2)$mcd.wt)))
## check for covMcd() consistency:
print(tablen <- table(vapply(L1, length, 1)))
plot(tablen)
print(iCommon <- Reduce(intersect, L1))
stopifnot(ii8 %in% iCommon)
##

