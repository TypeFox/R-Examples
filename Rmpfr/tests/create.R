require("Rmpfr")

### Simple basic examples of creation of   "mpfr"  objects

pi. <- Const("pi", prec = 260)
pi. # nicely prints 80 digits [260 * log10(2) ~= 78.3 ~ 80]

## These both failed (in mpfr2str(.)) with a seg.fault:
c(mpfr(1, prec=3), pi.)
m0 <- mpfr(numeric(), prec=64)
## print()ing / str() of 0-length mpfr
stopifnot(
    grepl("0 'mpfr' numbers", capture.output(    m0)),
    grepl("0 'mpfr' numbers", capture.output(str(m0))))


## This is TRUE for 0 and -0 :
Zero <- mpfr(c(0,1/-Inf), 20)
stopifnot(mpfrIs0(Zero), is.whole(Zero))
stopifnot(mpfr.is.0(Zero))# deprecated but must work
stopifnot(mpfr.is.integer(Zero))# deprecated but must work
Zero # the "-0" should print correctly
stopifnot(.getSign(Zero) == c(1,-1),
          sign(Zero) == 0,
	  identical(format(Zero, digits=1), c("0.", "-0.")))

## testing 'recycling'
b <- c(20,120,80, 60)
x <- mpfr(2^-(5:7), precBits = b);x

d.spec <- c(0,NA,NaN,Inf,-Inf)
(spec <- mpfr(d.spec, 3))
stopifnot(length(x) == 4, x[1] == x[4], getPrec(x) == b,
	  identical(is.na(spec), is.na(d.spec)),
	  identical(is.finite(spec), is.finite(d.spec)),
	  identical(is.infinite(spec), is.infinite(d.spec)),
	  identical(format(spec),
		    c("0.0", "NaN", "NaN", "Inf", "-Inf")),
	  mpfr(0.2, prec = 5:15, rnd.mode = "D") < 0.2)

x <- c(-12, 1:3 * pi)
sss <- mpfr(x, 100)
validObject(sss)
sss
sss2 <- sss * sss
stopifnot(identical(sss2, sss * x),
          identical(sss2, x * sss),
          sss ^ 2 == sss2)
## and go back {not sure if  identical() is guaranteed here, but it seems...}:
stopifnot(identical(x, as(sss, "numeric")))

(cs <- as(sss, "character"))

y <- c(0, 100,-10, 1.25, -2.5,
       x * c(1,100,1e5,1e20),
       x / 100^(1:4))
(Y <- mpfr(y, 100))
cbind(y, as.data.frame(.mpfr2str(Y, 20))[,c("exp","str")])

s <- mpfr(43208,  14)# low precision
eps8 <- 8 * .Machine$double.eps
## checking  mpfr -> character -> mpfr:
stopifnot(all.equal(y, as.numeric(format(Y, digits=20)), tol= eps8),
	  all.equal(Y, as(format(Y), "mpfr"), tol= eps8),
	  identical(sapply(1:5, formatMpfr, x=s),
		    c("4.e4", "4.3e4", "4.32e4", "43210.", "43208.")))


## More  character -> mpfr  checking :
## from   echo 'scale=200; 4*a(1)' | bc -l :
cpi <- "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196"
pi. <- Const("pi", prec=667)
stopifnot(cpi == format(mpfr(cpi, prec=667), digits=201),
          all.equal(pi., as(cpi, "mpfr")),
          all.equal(pi., as(cpi, "mpfr"), tol = 1e-200))

set.seed(17)
## Check double -> mpfr -> character -> double :
##  Unfortunately,  format(<mpfr>, .) -> .mpfr2str() triggers a memory bug
##  that I think is an MPFR library "mis-feature"
## 2011-02-09 -- bug *no longer* triggered !
rSign <- function(n) sample(c(-1,1), size = n, replace=TRUE)
N <- function(x) as.numeric(x)
ntry <- if(Sys.getenv("USER") == "maechler") 150 else 5
for(n in 1:ntry) {
    cat(if(n %% 10)"." else n)
    x. <- rSign(100) * rlnorm(100)
    prec <- rpois(1, 110); digs <- floor(0.95*(prec / log2(10)))
    X. <- mpfr(x., precBits = prec)
    stopifnot(all.equal(x., N(format(X., digits=digs)), tol = eps8)
              , all.equal(x., N(log(exp(X.))), tol = 32*eps8)
    )
}; cat("\n")

stopifnot(identical(mpfr.is.0(X.),# deprecated but must work
		    mpfrIs0  (X.)))
X. <- X.[!mpfrIs0(X.)]
stopifnot(all( X./X. == 1)) # TRUE

u <- mpfr(as.raw(0:100))
z <- mpfr(1:12, 200)
z[z > 100] <- 100 # nothing done  (but used to fail)
z[] <- 0
stopifnot(0:100 == u, is(z,"mpfr"), mpfrIs0(z),
	  all.equal(u, mpfr(0:100, prec = 8), tol = 0),
	  0:1 == mpfr(1:2 %% 2 == 0))
z[3] <- Const("pi",200)
## z has length 12 -- now extend it:
z[15:17] <- 1/mpfr(10:12, 100)
stopifnot(all.equal(z[1:4], c(0,0,pi,0), tol = 1e-15), validObject(z),
	  all.equal(z[13:17], c(NaN,NaN, 1/(10:12)), tol = 1e-15))

## These seg.faulted (each via different R -> C interface) in the past:
assertError <- tools::assertError
assertError( pp <- Const("pi", prec   = 1e11) )
assertError( mpfr("123.456",  precBits= 1e11) )
assertError( mpfr(as.bigz(3), precBits= 1e11) )
