svcm <- function(Y, X, vsize = c(1, 1, 1), knots = c(10, 10, 10),
                 deg = c(1, 1, 1), opd = c(1, 1, 1), search = TRUE,
                 lambda.init = rep(1e-3, 3), method = "grid", type = "SEQ",
                 ...) {

  t0 <- Sys.time()
  cl <- match.call()
  cat(paste("\nStarting ....\n", deparse(cl, width.cutoff=500),"\n"))
  ##avoid scoping conflicts
  svcm.env <- environment()
  environment(predictsvcm) <- svcm.env
  environment(targetsvcm) <- svcm.env
  environment(EDofTP) <- svcm.env

  dims <- dim(Y)
  r <- nrow(X)
  if (r != dims[length(dims)])
    stop("\nLast dimension of Y must match first dimension of X!\n")
  p <- ncol(X)
  n <- prod(dims)                       ##needed for GCV
  dims <- dims[-length(dims)]
  ndims <- length(dims)
  vsize <- vsize[1:ndims]
  
  X <- Matrix(X, sparse=TRUE)
  if (class(X) == "ddiMatrix")  X <- as(as(X, "dgeMatrix"), "dgCMatrix")
  Y <- aperm(Y, c(ndims + 1, 1:ndims))
  
  #####################################################
  ##     preparation of design and penalty matrices   #
  #####################################################

  ##compute absolute coordinates (voxel centers)
  x.coords <- (1:dims[1])*vsize[1] - vsize[1]/2
  y.coords <- (1:dims[2])*vsize[2] - vsize[2]/2
  
  ##evaluate basis functions
  B.x <- Matrix(splineDesign(calknots(dims[1], vsize[1], deg[1], knots[1]),
                             x.coords, ord = deg[1] + 1), sparse=TRUE)
  B.y <- Matrix(splineDesign(calknots(dims[2], vsize[2], deg[2], knots[2]),
                             y.coords, ord = deg[2] + 1), sparse=TRUE)
  
  ##number of coefficients in x- and y-direction, i.e. knots[i] + deg[i] - 1 
  r.x <- ncol(B.x)
  r.y <- ncol(B.y)
  
  ##compute difference matrices
  D.x <- Matrix(diff(diag(r.x), differences = opd[1]), sparse=TRUE)
  D.y <- Matrix(diff(diag(r.y), differences = opd[2]), sparse=TRUE)
  
  if (ndims == 3) {
    z.coords <- (1:dims[3])*vsize[3] - vsize[3]/2
    B.z <- Matrix(splineDesign(calknots(dims[3], vsize[3], deg[3], knots[3]),
                               z.coords, ord = deg[3] + 1) , sparse=TRUE)
    r.z <- ncol(B.z)
    D.z <- Matrix(diff(diag(r.z), opd[3]), sparse=TRUE)
  }
  
  ##adjust for 2d or 3d depending on SVCM type
  if (type == "TP") {
    
    if(ndims == 2) {    
      B <- as(B.y %x% B.x %x% X, "dgCMatrix")
      ##compute penalty matrices
      P.x <-  crossprod(as(Diagonal(r.y) %x% D.x %x% Diagonal(p),"dgCMatrix"))
      P.y <- crossprod(as(D.y %x% Diagonal(r.x*p), "dgCMatrix"))
      if (length(lambda.init) == 1) {
        P.whole <- as(P.x + P.y, "dsCMatrix")
        rm(P.x, P.y)
      }    
    } else if (ndims == 3) {
      B <- as(B.z %x% B.y %x% B.x %x% X,"dgCMatrix")
      ##compute penalty matrices
      P.x <- crossprod(as(Diagonal(r.z*r.y) %x% D.x %x% Diagonal(p),
                          "dgCMatrix"))
      P.y <- crossprod(as(Diagonal(r.z) %x% D.y %x% Diagonal(r.x*p),
                          "dgCMatrix"))
      P.z <- crossprod(as(D.z %x% Diagonal(r.y*r.x*p), "dgCMatrix"))
      if (length(lambda.init) == 1) {
        P.whole <- as(P.x + P.y + P.z, "dsCMatrix")
        rm(P.x, P.y, P.z)
        }    
    }
    ##matrices involved in the least squares estimation
    ##dsC2env requires U-form:
    BB <- crossprod(B)
    y <- as.vector(array(Y, c(r*dims[1], dims[2:ndims])))
    RHS <- crossprod(B, y)
    
  } else if (type == "SEQ") {
    
    ##compute penalty matrices
    P.x <- crossprod(D.x)
    P.y <- crossprod(D.y)
    ##matrices involved in the least squares estimation     
    LS.X <- solve(crossprod(X)) %*% as.matrix(t(X))
    BB.x <- crossprod(B.x)
    BB.y <- crossprod(B.y)
    
    if (ndims == 3) {    
      BB.z <- crossprod(B.z)
      P.z <- crossprod(D.z)
    }
  }
  ##clean up
  rm(x.coords, y.coords, D.x, D.y)
  if (ndims == 3) rm(z.coords, D.z)

  
  ##################################################
  ##    optimization of the smoothing parameter    #
  ##################################################
  
  if (search) {

    GCVtab <- NULL

    if (length(lambda.init) == 1) {     ## global       
      if (method == "grid") {
        opt <- cleversearch(targetsvcm, ...)
      } else {
        opt <- optimize(f = targetsvcm, ...)
        names(opt) <- c("par", "value")
      }
    } else {                            ##dimension-specific 
      if (method == "grid") {
        opt <- cleversearch(targetsvcm, ...)
      } else {
        opt <- optim(lambda.init, targetsvcm, method = method, ...)
      }
    }
    colnames(GCVtab) <- c(paste("lambda", 1:length(lambda.init), sep = ""),
                          "GCV")
    opt$GCVtab <- GCVtab
    predictsvcm(opt$par, type = type)
    t1 <- Sys.time()
    opt$time <- paste(t1 - t0, attributes(t0 - t1)$units)
    
  } else {
    
    predictsvcm(lambda.init, type = type)
    t1 <- Sys.time()
    opt <- list(time = paste(t1 - t0, attributes(t0 - t1)$units),
                par = lambda.init)
  }
  
  #########################################################################
  ##  compute return values and reorder to (nx,ny,nz,r) and (rx,ry,rz,p)  #
  ######################################################################### 

  if (type == "TP") {
    eta <- array(as.matrix(eta), c(r, dims))
    if (ndims == 2) {
      ##use tensorproduct structure: (B' %x% I_p)vec(A) = vec(AB)
      beta <- aperm(array(as.matrix(tcrossprod(Matrix(A, nrow = p),
                                               Matrix(B.y %x% B.x))),
                          c(p, dims)), c((1:ndims) + 1, 1))
      A <- array(A, c(p, r.x, r.y)) 
    } else if (ndims == 3) {
      beta <- aperm(array(as.matrix(tcrossprod(Matrix(A, nrow = p),
                                               Matrix(B.z %x% B.y %x% B.x))),
                          c(p, dims)), c((1:ndims) + 1, 1))
      A <- array(A, c(p, r.x, r.y, r.z))
    }
  } else if (type == "SEQ") {
    beta <- SEQpsi(diag(p), A, 1)
    beta <- SEQpsi(as.matrix(B.x), beta, 2)
    beta <- SEQpsi(as.matrix(B.y), beta, 3)
    if (ndims == 3) {
      beta <- SEQpsi(as.matrix(B.z), beta, 4)
    }
    beta <- aperm(beta, c((1:ndims) + 1, 1))
    Y <- aperm(Y, c((1:ndims) + 1, 1))
  }
  
  eta <- aperm(eta, c((1:ndims) + 1, 1))

  invisible(list(fitted = eta, effects = beta, coeff = A,
              knots = knots, deg = deg, opd = opd, vsize = vsize,
              type = type, call = cl, opt = opt))
  
}
