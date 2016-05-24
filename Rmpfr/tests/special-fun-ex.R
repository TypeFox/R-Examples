stopifnot(require("Rmpfr"))
(doExtras <- Rmpfr:::doExtras())

all.eq.finite <- function(x,y, ...) {
    ## x = 'target'   y = 'current'
    if(any(is.finite(y[!(fx <- is.finite(x))])))
	return("current has finite values where target has not")
    if(any(is.finite(x[!(fy <- is.finite(y))])))
	return("target has finite values where current has not")
    ## now they have finite values at the same locations
    all.equal(x[fx], y[fy], ...)
}
n <- 1000
head(x <- mpfr(0:n, 100) / n)

stopifnot(range(x) == 0:1
	  ,all.equal(as.numeric(j0(x)),
		     besselJ(as.numeric(x), 0), tol = 1e-14)
	  ,all.equal(as.numeric(j1(x)),
		     besselJ(as.numeric(x), 1), tol = 1e-14)
	  ,all.equal(as.numeric(y0(x)),
		     besselY(as.numeric(x), 0), tol = 1e-14)
	  ,all.equal(as.numeric(y1(x)),
		     besselY(as.numeric(x), 1), tol = 1e-14)
	  )

### pnorm() -> erf() :
u <- 7*x - 2
stopifnot(all.equal(pnorm(as.numeric(u)),
		    as.numeric(pnorm(u)), tol = 1e-14))
## systematic random input testing:
set.seed(101)
if(doExtras) {
    nSim <- 50
    n2 <- 100
} else {
    nSim <- 10
    n2 <- 64
}
for(n in 1:nSim) {
    N <- rpois(1, lambda=n2)
    N3 <- N %/% 3
    x <- c(rnorm(N-N3), 10*rt(N3, df=1.25))# <- some large values
    m <- rnorm(N, sd = 1/32)
    s <- rlnorm(N, sd = 1/8)
    cEps <- .Machine$double.eps
    for(LOG in c(TRUE,FALSE))
	for(L.T in c(TRUE,FALSE)) {
	    p. <- pnorm( x, m=m,sd=s, log.p=LOG, lower.tail=L.T)
	    stopifnot(all.equal(p., pnorm(mpfr(x, precBits= 48), m=m,sd=s, log.p=LOG, lower.tail=L.T),
				tol = 128 * cEps))
	    stopifnot(all.equal(p., pnorm(mpfr(x, precBits= 60), m=m,sd=s, log.p=LOG, lower.tail=L.T),
				tol = 2 * cEps))
	}
    cat(".")
};cat("\n")
proc.time()


### Riemann's Zeta function:

## -- integer arguments --
stopifnot(all(mpfrIs0(zeta(-2*(1:100)))))

k.neg <- 2*(-100:0) - 1
Z.neg <- zeta(k.neg)
plot(k.neg, abs(as.numeric(Z.neg)), type = "l", log="y")

Pi <- Const("pi", 128L)

## confirm published value of Euler's gamma to 100 digits
pub.g <-
    paste("0.5772156649", "0153286060", "6512090082", "4024310421", "5933593992",
	  "3598805767", "2348848677", "2677766467", "0936947063", "2917467495",
	  sep="")

## almost =
our.g <- Const("gamma", log2(10) * 100) # 100 digits
(ff.g <- .mpfr2str(our.g))


M <- function(x) mpfr(x, 128L)
stopifnot(all.equal(zeta( 0), -1/2,      tol = 2^-100)
	  , all.equal(zeta(-1), -1/M(12),  tol = 2^-100)
	  , all.equal(zeta(-3),  1/M(120), tol = 2^-100)
	  ## positive ones :
	  , all.equal(zeta(2),  Pi^2/6,   tol = 2^-100)
	  , all.equal(zeta(4),  Pi^4/90,  tol = 2^-100)
	  , all.equal(zeta(6),  Pi^6/945, tol = 2^-100)
	  )

### Exponential Integral Ei(.)
curve(Ei, 0,5, n=5001)
if(mpfrVersion() >= "3") { ## only available since MPFR 3.0.0
  ### Airy function Ai(.)
  curve(Ai, -10, 5, n=5001); abline(h=0,v=0, col="gray", lty=3)
}

### Utilities  hypot(), atan2() : --- TODO !

## beta(), lbeta()
## ---------------
## The simplistic "slow" versions:
B  <- function(a,b) { a <- as(a, "mpfr"); b <- as(b, "mpfr"); gamma(a)*gamma(b) / gamma(a+b) }
lB <- function(a,b) { a <- as(a, "mpfr"); b <- as(b, "mpfr"); lgamma(a)+lgamma(b) - lgamma(a+b) }

## For partly *integer* arguments
Bi1 <- function(a,b) 1/(a*chooseMpfr(a+b-1, a)) # a must be integer >= 0
Bi2 <- function(a,b) 1/(b*chooseMpfr(a+b-1, b)) # b must be integer >= 0

x <- 1:10 + 0 ; (b10 <- mpfr(x, 128L))

stopifnot(all.equal(	B(1,b10),  1/x),
	  all.equal(	B(2,b10),  1/(x*(x+1))),
	  all.equal( beta(1,b10),  1/x),
	  all.equal( beta(2,b10),  1/(x*(x+1))),
	  TRUE)

x <- -10:10 + 0; X <- mpfr(x, 128L)
stopifnot(Bi1(1,X) == (B1x <- Bi2(X,1)),
	  Bi1(2,X) == (B2x <- Bi2(X,2)),
	  Bi1(3,X) == (B3x <- Bi2(X,3)),
	  all.equal(B1x,  1/x,               tol= 4e-16) ,
	  all.equal(B2x,  1/(x*(x+1)),       tol= 8e-16) ,
	  all.equal(B3x,  2/(x*(x+1)*(x+2)), tol=16e-16) ,
	  ## these the "poles" are all odd i.e. result in { +Inf / -Inf / NaN}
	  ## are all "ok" {e.g. 1/(x*(x+1)) gives (-Inf, Inf) for x = -1:0 }
	  all.eq.finite(beta(1,X),  1/x) ,
	  all.eq.finite(beta(X,2),  1/(x*(x+1))) ,
	  all.eq.finite(beta(3,X),  2/(x*(x+1)*(x+2)), tol=16e-16) ,
	  TRUE)

## (a,b)  *both* integer, one negative:
for(i in (-20):(-1)) {
    cat(i,":\n")
    a <- mpfr(i, 99)
    i1 <- i+1
    b. <- seq_len(-i1)
    Bab <- beta(a, b.)
    stopifnot(is.nan(beta(a, (i1:0))), is.nan(lbeta(a, (i1:0))),
	      all.equal(Bab, Bi2(a, b.),             tol=1e-20),
	      all.equal(lbeta(a, b.), log(abs(Bab)), tol=1e-20))
}

## (a,b) all positive
c10 <- b10 + 0.25
for(a in c(0.1, 1, 1.5, 2, 20)) {
    stopifnot(all.equal( B(a,b10), (bb <-  beta(a, b10))),
	      all.equal(lB(a,b10), (lb <- lbeta(a, b10))),
	      all.equal(lb, log(bb)),
	      all.equal( B(a,c10), (bb <-  beta(a, c10))),
	      all.equal(lB(a,c10), (lb <- lbeta(a, c10))),
	      all.equal(lb, log(bb)),
	      TRUE)
}

## However, the speedup is *not* much (50%) when applied to vectors:
stopifnot(validObject(xx <- outer(b10, runif(20))),
	  dim(xx) == c(length(b10), 20),
	  validObject(vx <- as(xx, "mpfr")), class(vx) == "mpfr", is.null(dim(vx)))
C1 <- replicate(10, system.time(bb <<- beta(vx, vx+2)))
C2 <- replicate(10, system.time(b2 <<-    B(vx, vx+2)))
summary(1000*C1[1,]) ##  80.3 {cmath-5, 2009}
summary(1000*C2[1,]) ## 125.1 { " }
stopifnot(all.equal(bb, b2))
## and for a single number, the speedup is a factor 3:
x1 <- vx[1]; x2 <- x1+2
system.time(for(i in 1:100) bb <- beta(x1, x2))# .27
system.time(for(i in 1:100) b2 <-    B(x1, x2))# .83

## a+b is integer <= 0, but a and b are not integer:
a <- b <- .5 + -10:10
ab <- data.matrix(expand.grid(a=a, b=b, KEEP.OUT.ATTRS=FALSE))
ab <- mpfr(ab[rowSums(ab) <= 0, ], precBits = 128)
stopifnot( beta(ab[,"a"], ab[,"b"]) == 0,
	  lbeta(ab[,"a"], ab[,"b"]) == -Inf)
## was  NaN  in Rmpfr <= 0.5-2

stopifnot(all.equal(6 * beta(mpfr(1:3,99), -3.), c(-2,1,-2), tol=1e-20))
## add more checks, notably for b (> 0)  above and below the "large_b" in
## ../src/utils.c :
bb <- beta(mpfr(1:23, 128), -23)
stopifnot(all.equal(bb, Bi1(1:23, -23), tol=1e-7))
                                        # Bi1() does not get high prec for small b
## can be written via rationals:  N / D :
bn <- c(330, -360, 468, -728, 1365, -3120, 8840, -31824,
        151164, -1007760, 10581480, -232792560)
bn <- c(rev(bn[-1]), bn)
bd <- 24* as.bigz(2 * 3 * 5 * 7 * 11) * 13 * 17 * 19 * 23
stopifnot(all.equal(bb, as(bn/bd,"mpfr"), tol=0))

stopifnot(all.equal(6 * beta(mpfr(1:3,	99), -3.),
			     c(-2,1,-2),	    tol=1e-20),
	  all.equal(   lbeta(mpfr(1:3, 128), -3.),
		    log(mpfr(c( 2,1, 2), 128) / 6), tol=1e-20))

## add more checks, notably for b (> 0)  above and below the "large_b" in
## ../src/utils.c :
bb <- beta(mpfr(1:23, 128), -23)
stopifnot(all.equal(bb, Bi1(1:23, -23), tol=1e-7))
					# Bi1() does not get high prec for small b
## can be written via rationals:  N / D :
bn <- c(330, -360, 468, -728, 1365, -3120, 8840, -31824,
	151164, -1007760, 10581480, -232792560)
bn <- c(rev(bn[-1]), bn)
bd <- 24* as.bigz(2 * 3 * 5 * 7 * 11) * 13 * 17 * 19 * 23
stopifnot(all.equal(bb, as(bn/bd,"mpfr"), tol=0))

## 2) add check for 'b' >  maximal unsigned int {so C code uses different branch}
two <- mpfr(2, 128)
for(b in list(mpfr(9, 128), mpfr(5, 128)^10, two^25, two^26, two^100)) {
    a <- -(b+ (1:7))
    stopifnot(a+b == -(1:7), # just ensuring that there was no cancellation
	      is.finite( B <-  beta(a,b)), ## was NaN ..
	      is.finite(lB <- lbeta(a,b)), ## ditto
	      all.equal(log(abs(B)), lB),
	      TRUE)
}

ee <- c(10:145, 5*(30:59), 10*(30:39), 25*(16:30))
b <- mpfr(2, precBits = 10 + max(ee))^ee # enough precision {now "automatic"}
stopifnot((b+4)-b == 4, # <==> enough precision above
	  b == (b. <- as(as(b,"bigz"),"mpfr")))
(pp <- getPrec(b.))# shows why b. is not *identical* to b.
system.time(Bb <- beta(-b-4, b))# 0.334 sec
if(dev.interactive())
    plot(ee, asNumeric(log(Bb)), type="o",col=2)
lb <- asNumeric(log(Bb))
## using  coef(lm(lb ~ ee))
stopifnot(all.equal(lb, 3.175933 -3.46571851*ee, tol = 1e-5))# 4.254666 e-6


bb <- beta(           1:4,   mpfr(2,99))
stopifnot(identical(bb, beta(mpfr(2,99), 1:4)),
	  all.equal((2*bb)*cumsum(1:4), rep(1, 4), tol=1e-20),
	  getPrec(bb) == 128)


cat('Time elapsed: ', proc.time(),'\n') # "stats"
if(!interactive()) warnings()
