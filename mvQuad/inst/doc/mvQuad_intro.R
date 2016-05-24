## ------------------------------------------------------------------------
  library(mvQuad)

  # create grid
  nw <- createNIGrid(dim=2, type="GLe", level=6)
  
  # rescale grid for desired domain
  rescale(nw, domain = matrix(c(1, 1, 2, 2), ncol=2))

  # define the integrand
  myFun2d <- function(x){
    x[,1]*exp(x[,2])
  }

  # compute the approximated value of the integral
  A <- quadrature(myFun2d, grid = nw)

## ------------------------------------------------------------------------
# via an user-defined function
  myRule.fun <- function(l){
    n <- seq(1, 2*l-1, by=2)/ (l*2)
    w <- rep(1/(l), l)

    initial.domain <- matrix(c(0,1), ncol=2)
    return(list(n=as.matrix(n), w=as.matrix(w), features=list(initial.domain=initial.domain)))
  }

  nw.fun <- createNIGrid(d=1, type = "myRule.fun", level = 10)
  print(data.frame(nodes=getNodes(nw.fun), weights=getWeights(nw.fun)))

## ------------------------------------------------------------------------
# via a text-file
  myRule.txt <- readRule(file=system.file("extdata", "oNC0_rule.txt", package = "mvQuad"))
  nw.txt <- createNIGrid(d=1, type = myRule.txt, level = 10)
  print(data.frame(nodes=getNodes(nw.txt), weights=getWeights(nw.txt)))

## ----echo=FALSE----------------------------------------------------------
  txt <- readLines(system.file("extdata", "oNC0_rule.txt", package = "mvQuad"))  
  for (i in 1:9) {
    cat(i,"\t", txt[i], "\n")
    if (i==9) {
      cat("   ... \n")    
    }
  }
  
  

## ----fig.show="hold", fig.width=4.5--------------------------------------
  nw <- createNIGrid(dim=2, type="cNC1", level=5, ndConstruction = "product")
  plot(nw, main="Example: Product-Rule")

## ----fig.show="hold", fig.width=4.5--------------------------------------
  nw <- createNIGrid(dim=2, type="cNC1", level=5, ndConstruction = "sparse")
  plot(nw, main="Example: Combination Technique")

## ---- fig.height=6, fig.width=6, echo=FALSE, message=FALSE---------------
  nw <- createNIGrid(dim=2, type="GHe", level=3)

  C = matrix(c(1,0.6,0.6, 2),2)
  m = c(2, 1)

  dmvnorm <- function(x,m,C){
    C.inv <- solve(C)
    tmp <- apply(x, MARGIN = 1, 
                 FUN = function(x){
                   1/sqrt((2*pi)^2 * det(C)) * exp(-0.5 * t(x-m)%*%C.inv%*%(x-m)) 
                 })
    return(unlist(tmp))
  }
  xx <- expand.grid(seq(-5,5,length.out=100),seq(-5,5,length.out=100))
  yy <- dmvnorm(xx, m, C)
  yy <- matrix(yy,100)

  par(mfrow=c(2,2), mar=c(2,2,3,1), oma=c(0,0,0,0))
  plot(nw, main="initial grid", xlim=c(-6,6), ylim=c(-6,6), pch=8)
  contour(seq(-5,5,length.out=100),seq(-5,5,length.out=100),yy,asp=1, labels="", add = T, col="gray")

  rescale(nw, m = m, C = C, dec.type = 0)
  plot(nw, main="no dec.", xlim=c(-6,6), ylim=c(-6,6), pch=8)
  contour(seq(-5,5,length.out=100),seq(-5,5,length.out=100),yy,asp=1, labels="", add = T, col="gray")
  

  rescale(nw, m = m, C = C, dec.type = 1)
  plot(nw, main="spectral dec.", xlim=c(-6,6), ylim=c(-6,6), pch=8)
  contour(seq(-5,5,length.out=100),seq(-5,5,length.out=100),yy,asp=1, labels="", add = T, col="gray")

  rescale(nw, m = m, C = C, dec.type = 2)
  plot(nw, main="Cholesky dec.", xlim=c(-6,6), ylim=c(-6,6), pch=8)
  contour(seq(-5,5,length.out=100),seq(-5,5,length.out=100),yy,asp=1, labels="", add = T, col="gray")

## ---- fig.height=6, fig.width=6, eval=FALSE------------------------------
#    nw <- createNIGrid(dim=2, type="GHe", level=3)
#  
#    # no rescaling
#    plot(nw, main="initial grid", xlim=c(-6,6), ylim=c(-6,6), pch=8)
#  
#    rescale(nw, m = m, C = C, dec.type = 0)
#    plot(nw, main="no dec.", xlim=c(-6,6), ylim=c(-6,6), pch=8)
#  
#  
#    rescale(nw, m = m, C = C, dec.type = 1)
#    plot(nw, main="spectral dec.", xlim=c(-6,6), ylim=c(-6,6), pch=8)
#  
#    rescale(nw, m = m, C = C, dec.type = 2)
#    plot(nw, main="Cholesky dec.", xlim=c(-6,6), ylim=c(-6,6), pch=8)

## ------------------------------------------------------------------------
  nw <- createNIGrid(dim=2, type="GHe", level=6, ndConstruction = "sparse")

  plot(nw)
  print(nw)
  dim(nw)
  size(nw)
  

