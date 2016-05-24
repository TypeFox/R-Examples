# library(rgl)

# Morris's screening method (Morris 1992, Campolongo 2007)
#
# Gilles Pujol 2006-2008
# Modified by Frank Weber (2016): support model functions
# returning a matrix or a 3-dimensional array.
#
# Sub-files:
# * morris_oat.R
# * morris_simplex.R
# * morris_sfd.R


ind.rep <- function(i, p) {
# indices of the points of the ith trajectory in the DoE
  (1 : (p + 1)) + (i - 1) * (p + 1)
}


morris <- function(model = NULL, factors, r, design, binf = 0, bsup = 1, scale = TRUE, ...) {
  
  # argument checking: factor number and names
  if (is.character(factors)) {
    X.labels <- factors
    p <- length(X.labels)
  } else if (is.numeric(factors)) {
    p <- factors
    X.labels <- paste("X", 1 : p, sep="")
  } else {
    stop("invalid argument \'factors\', waiting for a scalar (number) or a character string vector (names)")
  }
  
  # argument checking: number of repetitions
  if (length(r) == 1) {
    r.max <- r
  } else {
    r.max <- r[2]
    r <- r[1]
  }
  
  # argument checking: design parameters
  if (! "type" %in% names(design)) {
    design$type <- "oat"
    warning("argument \'design$type\' not found, set at \'oat\'")
  }
  if (design$type == "oat") {
    # one-at-a-time design
    if (! "levels" %in% names(design)) {
      stop("argument \'design$levels\' not found")
    }
    nl <- design$levels
    if (length(nl) == 1) nl <- rep(nl, p)
    if ("grid.jump" %in% names(design)) {
      jump <- design$grid.jump
      if (round(jump, 0) != jump) stop("grid.jump must be integer")
      if (length(jump) == 1) jump <- rep(jump, p)
    } else {
      jump <- rep(1, p)
      warning("argument \'design$grid.jump\' not found, set at 1")
    }
  } else if (design$type == "simplex") {
    # simplex-based design
    if (! "scale.factor" %in% names(design)) {
      stop("argument \'design$scale.factor\' not found")
    }
    h <- design$scale.factor
  } else {
    stop("invalid argument design$type, waiting for \"oat\" or \"simplex\"")
  }
  
  # argument checking: domain boundaries
  if (length(binf) == 1) binf <- rep(binf, p)
  if (length(bsup) == 1) bsup <- rep(bsup, p)
  
  # generation of the initial design
  if (design$type == "oat") {
    X <- random.oat(p, r.max, binf, bsup, nl, jump)
  } else if (design$type == "simplex") {
    X <- random.simplexes(p, r.max, binf, bsup, h)
  }
  
  # duplicated repetitions are removed
  X.unique <- array(t(X), dim = c(p, p + 1, r.max))
  X.unique <- unique(X.unique, MARGIN = 3)
  X <- matrix(X.unique, ncol = p, byrow = TRUE)
  colnames(X) <- X.labels
  r.unique <- nrow(X) / (p + 1)
  if (r.unique < r.max) {
    warning(paste("keeping", r.unique, "repetitions out of", r.max))
  }
  r.max <- r.unique
  
  # optimization of the design
  if (r < r.max) {
    ind <- morris.maximin(X, r)
    X <- X[sapply(ind, function(i) ind.rep(i, p)),]
  }
  
  # object of class "morris"
  x <- list(model = model, factors = factors, r = r, design = design,
            binf = binf, bsup = bsup, scale = scale, X = X, call =
              match.call())
  class(x) <- "morris"
  
  # computing the response if the model is given
  if (!is.null(x$model)) {
    response(x, other_types_allowed = TRUE, ...)
    tell(x)
  }
  
  return(x)
}


tell.morris <- function(x, y = NULL, ...) {
  id <- deparse(substitute(x))
  
  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }
  
  X <- x$X
  y <- x$y
  
  if (x$scale) {
    #X <- scale(X)
    #y <- as.numeric(scale(y))
    Binf <- matrix(x$binf, nrow = nrow(X), ncol = length(x$binf), byrow = TRUE)
    Bsup <- matrix(x$bsup, nrow = nrow(X), ncol = length(x$bsup), byrow = TRUE)
    X <- (X - Binf) / (Bsup - Binf) 
  }
  
  if (x$design$type == "oat") {
    x$ee <- ee.oat(X, y)
  } else if (x$design$type == "simplex") {
    x$ee <- ee.simplex(X, y)
  }
  
  assign(id, x, parent.frame())
}


print.morris <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y) && class(x$y) == "numeric") {
    cat("\nModel runs:", length(x$y), "\n")
    mu <- apply(x$ee, 2, mean)
    mu.star <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma <- apply(x$ee, 2, sd)
    
    out <- data.frame(mu, mu.star, sigma)
    rownames(out) <- colnames(x$ee)
    print(out)
  } else if (!is.null(x$y) && class(x$y) == "matrix") {
    cat("\nModel runs:", nrow(x$y), "\n")
    mu <- apply(x$ee, 3, function(M){
      apply(M, 2, mean)
    })
    mu.star <- apply(abs(x$ee), 3, function(M){
      apply(M, 2, mean)
    })
    sigma <- apply(x$ee, 3, function(M){
      apply(M, 2, sd)
    })
    out <- list("mu" = mu, "mu.star" = mu.star, "sigma" = sigma)
    print.listof(out)
  } else if (!is.null(x$y) && class(x$y) == "array") {
    cat("\nModel runs:", dim(x$y)[1], "\n")
    mu <- lapply(1:dim(x$ee)[4], function(i){
      apply(x$ee[, , , i, drop = FALSE], 3, function(M){
        apply(M, 2, mean)
      })
    })
    mu.star <- lapply(1:dim(x$ee)[4], function(i){
      apply(abs(x$ee)[, , , i, drop = FALSE], 3, function(M){
        apply(M, 2, mean)
      })
    })
    sigma <- lapply(1:dim(x$ee)[4], function(i){
      apply(x$ee[, , , i, drop = FALSE], 3, function(M){
        apply(M, 2, sd)
      })
    })
    names(mu) <- names(mu.star) <- names(sigma) <- dimnames(x$ee)[[4]]
    cat("----------------\nmu:\n\n")
    print.listof(mu)
    cat("----------------\nmu.star:\n\n")
    print.listof(mu.star)
    cat("----------------\nsigma:\n\n")
    print.listof(sigma)
  } else {
    cat("\n(empty)\n")
  }
}


plot.morris <- function(x, identify = FALSE, atpen = FALSE,
                        y_col = NULL, y_dim3 = NULL, ...) {
  if (!is.null(x$ee)) {
    if(class(x$y) == "numeric"){
      mu.star <- apply(x$ee, 2, function(x) mean(abs(x)))
      sigma <- apply(x$ee, 2, sd)
    } else if(class(x$y) == "matrix"){
      if(is.null(y_col)) y_col <- 1
      if(!is.null(y_dim3)){
        warning("Argument \"y_dim3\" is ignored since the model output is ",
                "a matrix")
      }
      mu.star <- apply(x$ee[, , y_col, drop = FALSE], 2, 
                       function(x) mean(abs(x)))
      sigma <- apply(x$ee[, , y_col, drop = FALSE], 2, sd)
    } else if(class(x$y) == "array"){
      if(is.null(y_col)) y_col <- 1
      if(is.null(y_dim3)) y_dim3 <- 1
      mu.star <- apply(x$ee[, , y_col, y_dim3, drop = FALSE], 2, 
                       function(x) mean(abs(x)))
      sigma <- apply(x$ee[, , y_col, y_dim3, drop = FALSE], 2, sd)
    }
    
    plot(mu.star, sigma, pch = 20, xlab = expression(mu^"*"),
         ylab = expression(sigma), ...)
    
    if (identify) {
      identify(mu.star, sigma, labels = colnames(x$ee), atpen = atpen)
    } else {
      text(mu.star, sigma, labels = colnames(x$ee), pos = 4)
    }
  }
}


plot3d.morris <- function(x, alpha = c(0.2, 0), sphere.size = 1,
                          y_col = NULL, y_dim3 = NULL) {
  spheres.rad <- max((apply(x$ee,2,max) - apply(x$ee,2,min)) / 100) * sphere.size
  color = "grey"
  cone.nfaces = 100
  
  if (!is.null(x$ee)) {
    if(class(x$y) == "numeric"){
      mu <- apply(x$ee, 2, mean)
      mu.star <- apply(x$ee, 2, function(x) mean(abs(x)))
      sigma <- apply(x$ee, 2, sd)
    } else if(class(x$y) == "matrix"){
      if(is.null(y_col)) y_col <- 1
      if(!is.null(y_dim3)){
        warning("Argument \"y_dim3\" is ignored since the model output is ",
                "a matrix")
      }
      mu <- apply(x$ee[, , y_col, drop = FALSE], 2, mean)
      mu.star <- apply(x$ee[, , y_col, drop = FALSE], 2, 
                       function(x) mean(abs(x)))
      sigma <- apply(x$ee[, , y_col, drop = FALSE], 2, sd)
    } else if(class(x$y) == "array"){
      if(is.null(y_col)) y_col <- 1
      if(is.null(y_dim3)) y_dim3 <- 1
      mu <- apply(x$ee[, , y_col, y_dim3, drop = FALSE], 2, mean)
      mu.star <- apply(x$ee[, , y_col, y_dim3, drop = FALSE], 2, 
                       function(x) mean(abs(x)))
      sigma <- apply(x$ee[, , y_col, y_dim3, drop = FALSE], 2, sd)
    }
  }
  
  if (requireNamespace("rgl", quietly = TRUE)){ rgl::open3d()}
  
  xmax <- max(mu.star)
  zmax <- max(max(sigma), xmax)
  
  n <- 100
  theta <- seq(from = 0, to = pi, length.out = n + 1)
  x <- rep(c(0, xmax, xmax), n)
  y <- as.numeric(rbind(rep(0, n), - xmax * cos(theta[-(n + 1)]),
                        - xmax * cos(theta[-1])))
  z <- as.numeric(rbind(rep(0, n), xmax * sin(theta[-(n + 1)]),
                        xmax * sin(theta[-1])))
  if (requireNamespace("rgl", quietly = TRUE)){ 
    rgl::triangles3d(x, y, z, color = color, alpha = alpha[1])
    x <- rep(c(0, xmax, xmax, 0), 2)
    y <- c(0, -xmax, -xmax, 0, 0, xmax, xmax, 0)
    z <- c(0, 0, zmax, zmax, 0, 0, zmax, zmax)
    
    rgl::quads3d(x, y, z, color = color, alpha = alpha[2])
    rgl::plot3d(mu.star, mu, sigma, type = "s", radius = spheres.rad, add = TRUE)
    rgl::axes3d()
    rgl::title3d(xlab="mu*", ylab="mu", zlab="sigma")
  }
}
