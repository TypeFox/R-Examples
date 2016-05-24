### Short comparison of the various cubature methods

## 2D-function to integrate (here: isotropic zero-mean Gaussian density)
f <- function (s, sigma = 5) exp(-rowSums(s^2)/2/sigma^2) / (2*pi*sigma^2)

## simple polygonal integration domain
disc.owin <- spatstat::disc(radius=5, centre=c(3,2), npoly=8)

## plot image of the function and integration domain
plotpolyf(disc.owin, f, xlim=c(-8,8), ylim=c(-8,8))


### Two-dimensional midpoint rule

testmidpoint <- function (eps, main=paste("2D midpoint rule with eps =",eps))
{
    plotpolyf(disc.owin, f, xlim=c(-8,8), ylim=c(-8,8), use.lattice=FALSE)    
    ## add evaluation points to plot
    with(spatstat::as.mask(disc.owin, eps=eps),
         points(expand.grid(xcol, yrow), col=m, pch=20))
    polyCub.midpoint(disc.owin, f, eps=eps)
}
testmidpoint(5)
testmidpoint(3)
testmidpoint(0.5)
testmidpoint(0.2)


### Product Gauss cubature using an increasing number of nodes

for (nGQ in c(1:5,10,20,60)) {
    cat("nGQ =", sprintf("%2i",nGQ), ": ",
        format(polyCub.SV(disc.owin, f, nGQ=nGQ), digits=16),
        "\n")
}

## 'rotation' affects location of nodes
opar <- par(mfrow=c(1,2))
polyCub.SV(disc.owin, f, nGQ=2, rotation=FALSE, plot=TRUE)
polyCub.SV(disc.owin, f, nGQ=2, rotation=TRUE, plot=TRUE)
par(opar)


### Line integration along the boundary for isotropic functions

polyCub.iso(disc.owin, f, center=c(0,0))   # see ?polyCub.iso


### Quasi-exact cubature of the bivariate Gaussian density
### using gpclib::tristrip and mvtnorm::pmvnorm()

if (requireNamespace("mvtnorm") && gpclibPermit()) {
    print(polyCub.exact.Gauss(disc.owin, mean=c(0,0), Sigma=5^2*diag(2),
                              plot=TRUE), digits=16)
}
