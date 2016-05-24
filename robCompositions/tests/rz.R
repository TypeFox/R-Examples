require(robCompositions)
require(MASS)
require(robustbase)
crnorm <- function(n, mu, Sigma){ 
  constSum(data.frame(exp(mvrnorm(n, mu, Sigma))))
}
sigGen <- function(p, d){
  x <- diag(p)
  x[upper.tri(x)] <- x[lower.tri(x)] <- d
  x
}
set.seed(1234)
x <- crnorm(50, rep(10,10), Sigma=sigGen(10,0.9))
x <- x+rnorm(ncol(x)*nrow(x), 10, 1)
lim <- quantile(as.numeric(as.matrix(x)), 0.05)
x[x < lim] <- 0
w <- x==0
dl <- rep(lim, ncol(x))

res1 <- impRZilr(x, dl=dl, method="lm")
res2 <- impRZilr(x, dl=dl, method="MM")
res3 <- impRZilr(x, dl=dl, method="pls", nComp="boot", verbose=TRUE)
res4 <- impRZilr(x, dl=dl, method="pls", nComp=rep(5,ncol(x)), verbose=TRUE)

epsilon <- 55*.Machine$double.eps
if(any(!c(isTRUE(is.logical((x[1,3]/x[1,4] - res1$x[1,3]/res1$x[1,4]) < epsilon)),
isTRUE(is.logical((x[1,3]/x[1,4] - res2$x[1,3]/res2$x[1,4]) < epsilon)),
isTRUE(is.logical((x[1,3]/x[1,4] - res3$x[1,3]/res3$x[1,4]) < epsilon)),
isTRUE(is.logical((x[1,3]/x[1,4] - res4$x[1,3]/res4$x[1,4]) < epsilon))))){
  stop("ratios are not preserved.")
}

check1 <- function(x, w, dl){
  check <- logical(ncol(x))
  for(i in 1:ncol(x)){
    check[i] <- any(x[w[,i],i] > dl[i])
  }
  if(any(check)){
    res <- FALSE
    stop("values are over detection limit")
  } else {
    res <- TRUE
    message("imputed data passed check below detection limit")
  }
  invisible(res)
}

if(all(!c(check1(res1$x, w, dl=dl),
check1(res2$x, w, dl=dl),
check1(res3$x, w, dl=dl),
check1(res4$x, w, dl=dl)))){
  stop("one method imputed above detection limit")
}

## HD data:
set.seed(1234)
x <- crnorm(50, rep(10,60), Sigma=sigGen(60,0.9))
lim <- 0.01
x[x < lim] <- 0
w <- x==0
dl <- rep(lim, ncol(x))

resHD1 <- impRZilr(x, dl=dl, method="pls", nComp="boot", verbose=TRUE)
resHD2 <- impRZilr(x, dl=dl, method="pls", nComp=rep(5,ncol(x)), verbose=TRUE, maxit=2)


data(arcticLake)
x <- arcticLake
## generate rounded zeros artificially:
x[x[,1] < 10, 1] <- 0
xia <- impRZilr(x, dl=c(10,44,0), eps=0.01, method="MM")
xia$x


#data(expenditures)
#xOrig <- x <- expenditures

## DL as negative values
#x[x < 350] <- - 350

#impCoda2(x)

## DL given:
#impCoda2(xOrig, dl=c(0, 350, 350, 350, 350))


#imp <- impCoda(x, method='roundedZero')
#imp2 <- alrEM(x, pos=2, dl=rep(5,3))$xImp


