# ############################################################################
#  Demo for the mvQuad-package
#
# Filename: mvQuad_demo.R
# Author: Constantin Weiser (weiserc@hhu.de)
# Creation date: 2015-10-12
# ############################################################################
# Description:
#   demonstrates the features of the package and his typical usage
#


## standard usage
  # create grid
  nw <- createNIGrid(dim=2, type="GLe", level=6)
  # rescale grid for required support
  rescale(nw, domain = matrix(c(1, 1, 2, 2), ncol=2))
  # summaries for grid
  print(nw)
  # plot the grid
  plot(nw)

  myFun2d <- function(x){
    x[,1]*exp(x[,2])
  }

  # apply quadrature to a concrete function
  quadrature(myFun2d, grid = nw)


## illustration product-rule vs. combination technique

  # create 3d grid with product-rule
  nw.p <- createNIGrid(dim=3, type="cNC1", level=4, level.trans = T)
  rescale(nw.p, domain = matrix(c(1, 1, 1, 2, 2, 2), ncol=2))
  print(nw.p)
  plot(nw.p)

  # create 3d grid with combination technique
  nw.s <- createNIGrid(dim=3, type="cNC1", level=4, ndConstruction = "sparse", level.trans = T)
  rescale(nw.s, domain = matrix(c(1, 1, 1, 2, 2, 2), ncol=2))
  print(nw.s)
  plot(nw.s)

  myFun3d <- function(x){
    x[,1] * exp(x[,2]) * sqrt(x[,3])
  }

  #true value: 8.54017036969688...
  quadrature(myFun3d, grid = nw.p)
  quadrature(myFun3d, grid = nw.s)


## rescaling methods
  nw <- createNIGrid(dim=2, type="GHe", level=3)

  C = matrix(c(1,0.9,0.9,3),2)
  m = c(-.5, .3)

  require(mvtnorm)
  xx <- expand.grid(seq(-5,5,length.out=100),seq(-5,5,length.out=100))
  yy <- dmvnorm(xx, m, C)
  yy <- matrix(yy,100)

  Ex <- function(x, m, C){
    x[,1] * dmvnorm(x, m, C)
  }

  par(mfrow=c(1,3), mar=c(9,4,4,1))
  rescale(nw, m = m, C = C, dec.type = 0)
  plot(nw, main="no decomposition", xlim=c(-6,6), ylim=c(-6,6), pch=8)
  points(nw$nodes[,1], nw$nodes[,2], pch=16, col="darkgray")
  contour(seq(-5,5,length.out=100),seq(-5,5,length.out=100),yy,asp=1, labels="", add = T, col="gray")
  quadrature(Ex, nw, m=m, C=C)


  rescale(nw, m = m, C = C, dec.type = 1)
  plot(nw, main="spectral decomposition", xlim=c(-6,6), ylim=c(-6,6), pch=8)
  points(nw$nodes[,1], nw$nodes[,2], pch=16, col="darkgray")
  contour(seq(-5,5,length.out=100),seq(-5,5,length.out=100),yy,asp=1, labels="", add = T, col="gray")
  quadrature(Ex, nw, m=m, C=C)
  par(xpd=TRUE)
  legend(x=0, y=-9, legend = c("initial", "rescaled"), pch=c(16, 8), col=c("darkgray", 1),
         xpd=TRUE, horiz = TRUE, xjust = 0.5, cex = 1.3)


  rescale(nw, m = m, C = C, dec.type = 2)
  plot(nw, main="Cholesky decomposition", xlim=c(-6,6), ylim=c(-6,6), pch=8)
  points(nw$nodes[,1], nw$nodes[,2], pch=16, col="darkgray")

  contour(seq(-5,5,length.out=100),seq(-5,5,length.out=100),yy,asp=1, labels="", add = T, col="gray")
  quadrature(Ex, nw, m=m, C=C)


  par(mfrow=c(1,1))
  nw <- createNIGrid(dim=2, type="GLe", level=4)
  plot(nw, xlim=c(-1,1), ylim=c(-1,1), pch=16, col="darkgray")
  rescale(nw, domain = rbind(c(-1, 1), c(-1, 1)))
  par(new=T)
  plot(nw, xlim=c(-1,1), ylim=c(-1,1), pch=8)
  legend(-0, -2, legend = c("initial", "rescaled"), col=c("darkgray", 1), pch=c(16, 8),
         xpd=TRUE, horiz = TRUE, xjust = 0.5, cex = 1.3)

## Self-defined quadrature rule (here: midpoint-rule)

# via a function
  myRule.fun <- function(l){
    n <- seq(1, 2*l-1, by=2)/ (l*2)
    w <- rep(1/(l), l)

    initial.domain <- matrix(c(0,1), ncol=2)
    return(list(n=as.matrix(n), w=as.matrix(w), features=list(initial.domain=initial.domain)))
  }

  nw.fun <- createNIGrid(d=1, type = "myRule.fun", level = 10)
  print(data.frame(nodes=getNodes(nw.fun), weights=getWeights(nw.fun)))


# via a text-file
  myRule.txt <- readRule(file=system.file("extdata", "oNC0_rule.txt", package = "mvQuad"))
  nw.txt <- createNIGrid(d=1, type = myRule.txt, level = 10)
  print(data.frame(nodes=getNodes(nw.txt), weights=getWeights(nw.txt)))


# via hard-coded rules
  nw.hc <- createNIGrid(d=1, type = "oNC0", level = 10)
  print(data.frame(nodes=getNodes(nw.hc), weights=getWeights(nw.hc)))

