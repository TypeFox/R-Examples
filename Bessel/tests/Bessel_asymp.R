#### Bessel functions  for (very) large  x/z  and/or  nu

#### TODO: move some stuff from ./IJKY.R to here  <<<<<

### The argument range is from David Scott, for K_nu(.) case:

library(Bessel)
xs <- 100*(1:7)
nus <- 450 + 10*(0:9)
d.xn <- expand.grid(nu = nus, x = xs)

M <- with(d.xn,
          cbind(K = besselK(x,nu), K_exp = besselK(x,nu, expon.scaled = TRUE),
                K_nA.2 = besselK.nuAsym(x, nu, log = TRUE, k.max=2),
                K_nA.3 = besselK.nuAsym(x, nu, log = TRUE, k.max=3),
                K_nA.4 = besselK.nuAsym(x, nu, log = TRUE, k.max=4))
          )
## Transform into nicely labelled 3d array :
A <- M
datt <- attr(d.xn, "out.attrs")
dim(A)      <- c(datt$dim,               ncol(M))
dimnames(A) <- c(datt$dimnames, list(colnames(M)))
A

## Compare the different approximation levels  k.max
stopifnot(
          all.equal(M[,3], M[,4], tol=1e-12)# 2.826 e-13
          ,
          all.equal(M[,4], M[,5], tol=2e-15)# 5.357 e-16
          ,
          all.equal(M[,"K"], exp(M[,5]), tol= 1e-12)# on log.scale: 2e-16 !
          )
