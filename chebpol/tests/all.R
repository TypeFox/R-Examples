library(chebpol)

cat("Have FFTW:",havefftw(),'\n')
set.seed(42)

a <- array(rnorm(24),c(x=2,y=3,z=4))
chebcoef(a)

if(havefftw()) {
# A long one-dimensional
f <- function(x) ifelse(x==0,0,sin(1/x))
ch <- chebappxf(f,dims=50000)
cat('Long test:',ch(0.03),f(0.03),ch(0.031),f(0.031),'\n')
}
f <- function(x) exp(-sum(x^2))

dims <- c(8,7,6)
ch <- chebappxf(f,dims)
s <- runif(3,-1,1)
ch(s)-f(s)
iv <- list(c(1,2),c(1,4),c(1,3))
ch <- chebappxf(f,dims,iv)
s <- c(1.4,2.3,1.9)
ch(s) - f(s)
a <- c(1.01,1.01,1.01) ; ch(a)- f(a)

sum(evalongrid(f,dims))
# vector valued function
g <- function(x) c(sum(x),prod(x),exp(-sum(x^2)))
gv <- evalongrid(g,dims)
sum(gv)
dim(gv)


chebknots(17)
chebknots(c(x=3,y=4,z=5))
# test chebappxg
## evenly spaced grid-points
su <- seq(0,1,length.out=10)
## irregularly spaced grid-points
s <- su^3
## create approximation on the irregularly spaced grid
ch <- Vectorize(chebappxg(exp(s),list(s)))
## test it:
r <- runif(1); cat('true:',exp(r),'appx:',ch(r),'\n')

#multivariate chebappxg
su <- seq(0,1,length.out=11)
grid <- list(su,su^2,su^3)
dims <- lapply(grid,length)

fv <- structure(apply(expand.grid(grid),1,f),dim=lapply(grid,length))
ch <- chebappxg(fv,grid)
s <- runif(3)
cat('true:',f(s),'appx:',ch(s),'\n')

# multi linear
s <- runif(3)
lip <- mlappx(fv,grid)
cat('true',f(s), 'appx:', ch(s), 'lip:',lip(s),'\n')

# test dct transform
a <- array(rnorm(24),c(2,3,4))
chebcoef(a,TRUE)

# uniform grid stuff
# Runge function
f <- function(x) 1/(1+25*x^2)
grid <- seq(-1,1,length.out=15)
val <- f(grid)
uc <- Vectorize(ucappx(val))
# and the Chebyshev
ch <- Vectorize(chebappxf(f,15))
# test it at 10 random points
t(replicate(10,{a<-runif(1,-1,1); c(arg=a, uc=uc(a), true=f(a), cheb=ch(a))}))

uc <- Vectorize(ucappxf(f,15,intervals=c(-1,1)))
a <- runif(1,-1,1)
uc(a); ch(a); f(a)

#polyharmonic splines
f <- function(x) 10/(10+sum(sqrt(x)))
knots <- matrix(runif(6000), 6)
phs <- polyh(f, knots, 3)
# test it in a random point
a <- runif(6)
f(a); phs(a)
phs <- polyh(f,knots,5)
phs(a)
phs <- polyh(f,knots,-20)
phs(a)
