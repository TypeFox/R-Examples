resolution <- function(X, svcmlist, fac) {
  
  t0 <- Sys.time()  
  dims <- dim(svcmlist$fitted)[-length(dim(svcmlist$fitted))]
  ndims <- length(dims)
  if (length(fac) != ndims) {
    stop("'fac' must be of same dimensionality as the coefficient surfaces.")
  }
  newdims <- as.integer(dims * fac)
  cat(paste("\nScaled dimensions will be", paste(newdims, collapse=" x "),
            ".\n"))
  newvsize <- svcmlist$vsize / fac
  r <- dim(X)[1]
  p <- dim(X)[2]
  A <- Matrix(svcmlist$coeff, nrow = p)
  
  ##compute knots used to estimate the SVCM 
  xknots <- calknots(dims[1], svcmlist$vsize[1], svcmlist$deg[1],
                     svcmlist$knots[1])
  yknots <- calknots(dims[2], svcmlist$vsize[2], svcmlist$deg[2],
                     svcmlist$knots[2])
  ##compute absolute coordinates of refined grid (voxel centers)
  nxcoords <- 0.5:(newdims[1] - 0.5) * newvsize[1]
  nycoords <- 0.5:(newdims[2] - 0.5) * newvsize[2]
  
  ##construction of the design matrix
  tB.x <- t(Matrix(splineDesign(xknots, nxcoords, ord = svcmlist$deg[1] + 1)))
  tB.y <- t(Matrix(splineDesign(yknots, nycoords, ord = svcmlist$deg[2] + 1)))
  if (ndims == 3) {
    zknots <- calknots(dims[3], svcmlist$vsize[3], svcmlist$deg[3],
                       svcmlist$knots[3])
    nzcoords <- 0.5:(newdims[3] - 0.5) * newvsize[3]
    tB.z <- t(Matrix(splineDesign(zknots, nzcoords,ord = svcmlist$deg[3] + 1)))
  }
  
  ##compute predicted values and effects using tensorproduct structures:
  ##(B' %x% X)vec(A) = vec(XAB) and (B' %x% I_p)vec(A) = vec(AB)
  if (ndims == 2) {
    eta <- aperm(array(as.matrix(X %*% A %*% (tB.y %x% tB.x)), c(r, newdims)),
                 c((1:ndims) + 1, 1))
    beta <- aperm(array(as.matrix(A %*% (tB.y %x% tB.x)), c(p, newdims)),
                  c((1:ndims) + 1, 1))
  } else if (ndims == 3) {
    eta <- aperm(array(as.matrix(X %*% A %*% (tB.z %x% tB.y %x% tB.x)),
                       c(r, newdims)),  c((1:ndims) + 1, 1))
    beta <- aperm(array(as.matrix(A %*% (tB.z %x% tB.y %x% tB.x)),
                        c(p, newdims)),  c((1:ndims) + 1, 1))
  }    
  t1 <- Sys.time()
  cat(paste("\nResolution scaling needed", t1 - t0, attributes(t0 - t1)$units,
            "\n"))

  invisible(list(fitted = eta, effects = beta))
  
}

