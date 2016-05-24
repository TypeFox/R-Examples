stopifnot(require("Rmpfr"))

x <- mpfr(0:7, 64)/7
mx <- x
dim(mx) <- c(4,2)
(m. <- mx) # "print"
m.[,2] <- Const("pi", 80)
m.[,] <- exp(mpfr(1, 90))
stopifnot(is(mx, "mpfrMatrix"), dim(mx) == c(4,2),
          is(m., "mpfrMatrix"), dim(m.) == dim(mx),
	  dim(is.finite(mx)) == dim(mx),
	  dim(is.nan(mx)) == dim(mx),
          getPrec(m.) == 90)

xx <- (0:7)/7
m.x <- matrix(xx, 4,2)
m2 <- mpfr(xx, 64); dim(m2) <- dim(m.x)
##
u <- 10*(1:4)
y <- 7 * mpfr(1:12, 80)
my <- y
dim(my) <- 3:4
m.y <- asNumeric(my)
stopifnot(all.equal(m2, mpfr(m.x, 64), tol=0), # not identical(..)
	  my[2,2] == 35,
	  my[,1] == 7*(1:3))

.N <- function(x) { if(!is.null(dim(x))) as(x,"array") else as(x,"numeric") }
noDN <- function(.) { dimnames(.) <- NULL ; . }
allEQ <- function(x,y) all.equal(x,y, tol=1e-15)

stopifnot(allEQ(m.x, noDN(.N(mx))),
	  allEQ(m.y, noDN(.N(my))),
	  allEQ(noDN(.N(my %*% mx)), m.y %*% m.x),
	  allEQ(noDN(.N(crossprod(mx, t(my)))), crossprod(m.x, t(m.y))),
	  allEQ(noDN(.N(tcrossprod(my, t(mx)))),
			tcrossprod(m.y, t(m.x))),
	  ##
	  identical(mx, t(t(mx))),
	  identical(my, t(t(my))),
	  ## matrix o vector .. even  vector o vector
	  identical(noDN(.N(my %*% 1:4)), m.y %*% 1:4 ),
	  identical(noDN(.N(my %*% my[2,])), m.y %*% .N(my[2,])),
	  identical( crossprod(1:3, my), 1:3 %*%   my),
	  identical(tcrossprod(1:4, my), 1:4 %*% t(my)),
	  identical(crossprod(y), t(y) %*% y),
	  identical(tcrossprod(y), y %*% t(y)),
	  identical(noDN(.N( crossprod(y))), crossprod(7 * 1:12)),
	  identical(noDN(.N(tcrossprod(y))),tcrossprod(7 * 1:12)),
	  identical(tcrossprod(1:3, u), noDN(.N(tcrossprod(1:3, as(u,"mpfr")))))
	  )

mx[3,1] <- Const("pi", 64)
stopifnot(allEQ(sum(mx[,1]), pi + 4/7))
m2 <- mx[c(1,4),]
stopifnot(dim(m2) == c(2,2), sum(m2) == 2)

## "mpfrArray" o "mpfr" :
Tmx <- array(TRUE, dim(mx), dimnames=dimnames(mx))
stopifnot(identical(Tmx, mx == (mx - mpfr(0, 10))),
	  identical(Tmx, mx - mpfr(1, 10) * mx == 0))
## subassignment, many kinds
mx[5] <- pi
mx[6] <- Const("pi",100)
stopifnot(validObject(mx), allEQ(mx[5], mx[6]),
	  getPrec(mx) == c(rep(64,5), 100, 64,64))

## %*% with vectors on LHS, ...
y <- t(2:4) # 1 x 3 matrix
m1 <-     (0:10)     %*% y
m2 <- mpfr(0:10, 50) %*% y
stopifnot((d <- m1 - m2) == 0, identical(dim(m1), dim(d)),
          m2 == m1, m1 == m2)

r <- 10*(0:4)
y <- t(2:6)
m1 <- 1:3 %*% y  %*% r
y. <- t(mpfr(2:6, 20))
m2 <- 1:3 %*% y. %*% r
stopifnot(m1 == m2, m1 - m2 == 0, identical(dim(m1), dim(m2)))

### Array (non-matrix) ---- indexing & sub-assignment :
A <- mpfrArray(1:24, prec = 96, dim = 2:4)
a <-     array(1:24,            dim = 2:4)
a.1 <- as(A[,,1], "array")
a1. <- as(A[1,,], "array")
A1. <- as(A[1,,], "mpfr")

stopifnot(all.equal(noDN(a.1), a[,,1], tol=0),
	  identical(A1., as.vector(A[1,,])),
	  ## arithmetic, subsetting etc:
	  allEQ(noDN(.N(A / A1.)), a/c(a1.)),
	  allEQ(noDN(.N(a / A1.)), a/c(a1.)),
          identical(noDN(A == 23), a == 23),
          identical(noDN(10 >= A), 10 >= a),
          identical(noDN(A <=  2), a <=  2),
          identical(noDN(A < 2.5), a < 2.5),
          identical(noDN(A !=  5), a !=  5),
          identical(A != 3, !(3 == A)),
          identical(-1 > A, A == 100),
          identical(noDN(A <= 0), a == pi)
	  )

A[1,2,3] <- Const("pi")
A[1, ,2] <- 1 / A[1,,2]
A

## check that A is "==" a  where
a <- array(1:24, 2:4); a[1,2,3] <- pi; a[1,,2] <- 1/a[1,,2]
stopifnot(allEQ(noDN(.N(A)), a),
          ## check aperm() methods :
          allEQ(noDN(.N(aperm(A))), aperm(a)),
          {p <- c(3,1:2); allEQ(noDN(.N(aperm(A,p))), aperm(a,p))},
          {p <- c(2:1,3); allEQ(noDN(.N(aperm(A,p))), aperm(a,p))})

## cbind() / rbind():
options(warn = 2)## no warnings here - ("exact recycling"):
validObject(m0 <- cbind(pi=pi, i = 1:6))
validObject(m1 <- cbind(a=Const("pi",60),i = 1:6, "1/mp" = 1/mpfr(1:3,70)))
validObject(m2 <- cbind(pi=pi, i = 1:2, 1/mpfr(1:6,70)))
validObject(n2 <- rbind(pi=pi, i = 1:2, 1/mpfr(1:6,70)))
stopifnot(is(m0,"matrix"), is(m1, "mpfrMatrix"), is(m2, "mpfrMatrix"),
          dim(m0) == c(6,2), dim(m1) == c(6,3), dim(m2) == c(6,3))
options(warn = 1)
suppressWarnings(eval(ex <- quote(m3 <- cbind(I=10, 1:3, inv=1/mpfr(2:3,80)))))
validObject(suppressWarnings(     n3 <- rbind(I=10, 1:3, inv=1/mpfr(2:3,80))))
stopifnot(identical(t(n2), m2),
          identical(t(n3), m3), validObject(m3),
          is(tryCatch(eval(ex), warning=function(.).), "warning"),
	  identical(cbind("A", "c"), matrix(c("A", "c"), 1,2)),
	  identical(rbind("A", 2),   matrix(c("A", "2"), 2,1)) )

## matrix(<mpfr>) works since 2015-02-28:
x <- mpfr(pi,64)*mpfr(2,64)^(2^(0:19))
(mx <- matrix(x, 4,5))

stopifnot(is(mx, "mpfrMatrix"),
    all.equal(matrix(0:19, 4,5),
              asNumeric(log2(log2(mx) - log2(Const("pi")))),
              tol = 1e-15)) # 64b-lnx: see 8.1e-17


cat('Time elapsed: ', proc.time(),'\n') # "stats"
if(!interactive()) warnings()
