###############################################################################
## misc helper functions
###############################################################################

#return a list of the interactionx and all its main effects
#for smooth terms, also add linear terms and interactions
allInvolvedTerms <- function(label, m) {
  trms <- strsplit(label, ":")[[1]]
  # add linear components to smooth terms
  if(any(grepl("sm(", trms, fixed = T))) {
    whichSm <- grep("sm(", trms, fixed = T)
    linTrms <- sapply(trms[whichSm], function(x) sub("sm(","lin(", x, fixed = T))
    # add interactions
    newInts <- expand.grid(c(trms, linTrms), c(trms, linTrms))
    newInts <- unique(apply(newInts, 1, paste, collapse =":"))
    trms <- unique(c(linTrms, trms, newInts))
  }
  return(unique(c(trms, label)[trms %in% names(m$predvars)]))
}


# center a basis matrix B s.t. its column space does not
# include polynomials of a vector x up to degree 'degree', or, if x is a matrix,
# s.t. B's new column space is orthogonal to x's column space
centerBase <- function(B, x, degree = 0, method = c("P", "QR")) {
  method <- match.arg(method)
  n <- nrow(B)
  if (degree == 0 && method == "QR") {
    # taken from mgcv/smooth.r l.1555 f.
    # see also S. Wood GAM:Introduction with R (ch. 1.8.1)
    # !this removes one column from B! (?!?)
    qrC <- qr(t(matrix(colSums(B), 1, ncol(B))))
    return(t(qr.qy(qrC, t(B)))[, -1, drop = FALSE])
  }
  else {
    # B_centered[, i] = resid(lm(B[, i] ~ 1 + x + ... + x^degree)) bzw
    #                  resid(lm(B[, i] ~ -1 + Xx))
    # bzw B_centered = (I - Xx (Xx'Xx)^-1 Xx') B
    Xx <- if(!is.matrix(x)) {
      tmp <- cbind(rep(1, n))
      if (degree > 0)
        tmp <- cbind(tmp, poly(x, degree))
      tmp
    } else x
    qrX <- qr(Xx)
    Bc <- qr.resid(qrX, B)
    return(Bc)
  }
}


# assign ... and only ... to the environment of function f
customFunctionEnv <- function(f, ...) {
  environment(f) <- new.env(parent = globalenv())
  dots <- list(...)
  nm <- names(dots)
  sapply(nm, function(x) {
    assign(x, dots[[x]], envir = environment(f))
  })

  return(f)
}

# extract models visited by the SSVS
getModels <- function(x, thresh = 0.5) {
  models <- unlist(lapply(x$samples$pV1,
    function(x) {
      apply(1 *(x > thresh), 1, paste, collapse ="")
    }))
  t <- table(models)
  c(sort(t/sum(t), decreasing = T))
}

#check whether points are inside 2-d polygon defined by vertices
#(assumes non-closed polygon, i.e head(vertices, 1)!= tail(vertices, 1))
insidePoly <- function(points, vertices) {

  #angle between two vectors, enforce (-180 deg, 180 deg)
  angle <- function(v1, v2) {
    ang <- (atan2(v1[2], v1[1]) - atan2(v2[2], v2[1]))
    while(ang > pi) ang <- ang - 2 * pi
    while(ang < -pi) ang <- ang + 2 * pi
    return(ang)
  }

  # idea: for a point inside the polygon, the sum of angles between
  # x and adjacent vertices has to be 2\pi (and 0 if x is outside ?)
  getAngleSum <- function(p, vertices) {
    vertexPairs <- cbind(1:nrow(vertices), c(2:nrow(vertices), 1))
    #check: p on a vertex?
    if(paste(p, collapse =";") %in% apply(vertices, 1, paste, collapse =";")) {
      return(2 * pi)
    }

    a <- apply(vertexPairs, 1, function(ind) {
      e1 <- vertices[ind[1],]- p
      e2 <- vertices[ind[2],] - p
      return(angle(as.numeric(e1), as.numeric(e2)))
    })
    #if point is exactly on hull, one angle must be |180|deg
    if(any(sapply(abs(a), identical, pi))) return(2 * pi)

    return(sum(a))
  }
  res <- apply(as.matrix(points), 1, function(p) {
    #all.equal(getAngleSum(p, vertices), 2 * pi)
    getAngleSum(p, as.matrix(vertices)) > pi
  })
  return(res)
}

# find starting values via QR-based IWLS
#' @import stats
iwls.start <- function(X, y, family, weights = rep(1, n),
  offset = rep(0, n), steps = 5) {
  n <- length(y)
  p <- NCOL(X)
  if(family == 2) {
    w <- sqrt(var <- mu <- drop(y + 0.1))
    eta <- log(mu)
    varBeta <- 4
  }
  if(family == 1) {
    mu <- drop((weights * y + 0.5)/(weights + 1))
    w <- sqrt(var <- weights * mu *(1-mu))
    eta <- qlogis(mu)
    varBeta <- 2
  }

  iwls.step <- function(coef, yWork, XWork) {
    coef <- qr.coef(qr(XWork), yWork)
    eta <- X %*% coef

    if(family == 2) {
      w <- sqrt(var <- mu <- drop(exp(eta + offset)))
    }
    if(family == 1) {
      mu <- drop(plogis(eta + offset))
      w <- sqrt(var <- weights * mu *(1-mu))
    }

    yWork <<- c(w *(eta + (y-mu)/var), rep(0, p))
    XWork <<- rbind(w * X, diag(1/sqrt(varBeta), p))

    return(coef)
  }

  yWork <- c(w *(eta + (y-mu)/var), rep(0, p))
  XWork <- rbind(w * X, diag(1/sqrt(varBeta), p))
  coef <- rep(0, p)
  for(i in 1:steps) (coef <- iwls.step(coef, yWork, XWork))

  return(coef)
}


# Return basis for penalized part of a function
# based on a matrix of basis functions B and associated penalty K
mmDesign <- function(B, K, tol = 1e-10, rankZ =.995) {
  ek <- eigen(K, symmetric = TRUE)
  nullvals <- ek$values < tol
  colsZ <- if(rankZ < 1) {
    #use at least 3 but at most so many eigenvectors that rankZ * 100% of info
    #in K is represented
    max(3, min(ncol(B), min(which( cumsum(ek$values[!nullvals]) /
        sum(ek$values[!nullvals]) >= rankZ ))))
  } else {
    rankZ
  }
  P <- t(1/sqrt(ek$values[1:colsZ]) * t(ek$vectors[, 1:colsZ]))
  return(B %*% P)
}

# Return orthogonalized low-rank centered basis for penalized part of a function
# based on a matrix of basis functions B with associated covariance Cov or
# on the implied (nonstationary) covariance of the function C = B Cov B'
orthoDesign <- function(B = NULL, Cov = diag(NCOL(B)),
  C = B %*%(Cov %*% t(B)), tol = 1e-10, rankZ =.995, rank = NULL) {
  ## do a full eigen decomposition only if there's no info on the rank of C
  ## i.e if neither B nor rank are supplied, otherwise use a truncated svd
  ## via lanczos-iteration (10-50 times faster)
  eC <- if(all(is.null(B), is.null(rank))) {
    eigen(C, symmetric = TRUE)
  } else {
    if(!is.null(B)) rank <- min(NROW(B), NCOL(B))
    if(is.null(rank)) rank<-qr(C)$rank
    try(irlba(C, rank, rank))
  }
  if(class(eC) == "try-error") {
    sv <- svd(C, rank, 0)
    eC <- list(values = sv$d, vectors = sv$u)
  }

  nullvals <- eC$values < tol
  colsZ <- if(rankZ < 1) {
    max(3, min(which( cumsum(eC$values[!nullvals])/sum(eC$values[!nullvals]) >
        rankZ )))
  } else {
    rankZ
  }
  colsZ <- min(colsZ, sum(!nullvals))
  use <- 1:colsZ
  return(t(sqrt(eC$values[use])* t(eC$vectors[, use])))
}

# make a P-spline basis with K basis functions
psBasis <- function(x, K = length(x), spline.degree = 3, diff.ord = 2,
  knots = NULL) {
  if(is.null(knots)) {
    knots.no <- K - spline.degree + 1
    xl <- min(x)
    xr <- max(x)
    xmin <- xl-(xr-xl)/100
    xmax <- xr +(xr-xl)/100
    dx <- (xmax-xmin)/(knots.no-1)
    knots <- seq(xmin-spline.degree * dx, xmax + spline.degree * dx, by = dx)
  }
  X <- splines::spline.des(knots, x, spline.degree + 1, outer.ok = TRUE)$design
  P <- diag(K) #precision
  if(diff.ord>0) {
    for(d in 1:diff.ord) P <- diff(P)
    P <-crossprod(P)
  }
  return(list(X = X, P = P, knots = knots, K = K, spline.degree = spline.degree,
    diff.ord = diff.ord))
}

#apply aggregate-function and quantiles to f
summarizeF <- function(f, aggregate, quantiles) {
  agg <- 	if(is.null(aggregate)) {
    NULL
  } else {
    apply(f, 1, aggregate)
  }
  quant <- if(is.null(quantiles)) {
    NULL
  } else {
    drop(t(apply(f, 1, quantile, probs = quantiles, na.rm = TRUE)))
  }
  res <- cbind(agg, quant)
  nm <- list(
    a = if(!is.null(aggregate)) {
      if(NCOL(agg)== 1) {
        "eta"
      } else {
        paste("eta.", 1:NCOL(agg), sep ="")
      }
    } else NULL ,
    qu = if(!is.null(quantiles)) {
      paste(quantiles * 100, "percentile", sep ="")
    } else NULL)
  colnames(res) <- unlist(nm)
  return(res)
}

# make a 2D thin plate spline base with K basis functions
# nd: boolean vector of NonDuplicated locations
# K: number of basis functions/knots
# knots: a matrix of coordinates for the knots
tp2DBasis <- function(coords, K = NULL, nd = rep(TRUE, length(coords)),
  knots = NULL) {
  get2dKnts <- function(coords, K) {
    # get good knot locations from 5 repeated runs of
    # a space-filling algorithm
    # TODO: use sampled points if coords is large?
    cl <-lapply(1:5, function(i) {
      c <- clara(coords, k = K, rngR = TRUE)
      return(list(obj = c$objective, kn = c$medoids))
    })
    return(cl[[which.min(sapply(cl, "[[", "obj"))]]$kn)
  }

  if(is.null(knots) & is.null(K)) stop("need to specify either K or knots.")
  if(is.null(knots)) knots <- get2dKnts(coords[nd,, drop = F], K)
  if(!is.null(knots)) K <- nrow(knots)

  rK <- dist(knots)
  P <-  as.matrix(rK^2 * log(rK))

  r1 <- outer(coords[, 1], knots[, 1], "-")
  r2 <- outer(coords[, 2], knots[, 2],"-")
  rX <- sqrt(r1^2 + r2^2)
  X <- rX^2 * log(rX)
  X[rX == 0] <- 0
  return(list(X = X, P = P, knots = knots, K = K))
}


#get symmetric A^-1/2
matrixInvRoot <- function(A, tol = 1e-10) {
  sva <- svd(A)
  Positive <- sva$d > max(tol * sva$d[1L], 0)
  if (all(Positive)) {
    return(sva$v %*% (1/sqrt(sva$d) * t(sva$u)))
  }
  else{
    return(sva$v[, Positive, drop = FALSE] %*% ((1/sqrt(sva$d[Positive]) *
        t(sva$u[, Positive, drop = FALSE]))))
  }
}


#compute sqrt(trace(x'x)/nrow(x))
frob <- function(x) {
  trace <- sum(sapply(1:NCOL(x),
    function(i) crossprod(x[, i])))
  return(sqrt(trace/nrow(x)))
}
#scale x s.t. it has ||x||_F = 1/factor
scaleMat <- function(x, factor = 2) {
  return(x/(factor * frob(x)))
}



safeDeparse <- function(expr) {
  ret <- paste(deparse(expr), collapse ="")
  #rm whitespace
  gsub("[[:space:]][[:space:]]+", " ", ret)
}

## Fast partial SVD by implicitly-restarted Lanczos bidiagonalization
## Arguments:
#@A  A double-precision real or complex matrix or real sparse matrix
#@nu  Number of desired left singular vectors nv 	Number of desired right
# singular vectors
#@adjust  Number of extra approximate singular values to compute to enhance
#  convergence
#@aug  "ritz" for Ritz "ham" for harmonic Ritz vector augmentation
#@sigma  "ls" for largest few singular values, "ss" for smallest
#@maxit  Maximum number of iterations
#@m_b  Size of the projected bidiagonal matrix
#@reorth Either 1 or 2: full (2) or one-sided (1) reorthogonalization
#@tol  Convergence is determined when ||A * V - U * S|| <= tol *||A||,
#  where ||A|| is approximated by the largest singular value of all projection
#  matrices.
#@V  Optional matrix of approximate right singular vectors
#
## Value:
#@d min (nu, nv) approximate singular values
#@u nu approximate left singular vectors
#@v nv approximate right singular vectors
# authors:  Jim Baglama and Lothar Reichel
# code taken from http://www.rforge.net/irlba/ (version 0.1.1)
irlba <- function (A, nu = 5, nv = 5, adjust = 3, aug ="ritz", sigma ="ls",
  maxit = 1000, m_b = 20, reorth = 2, tol = 1e-6, V = NULL) {

  # ---------------------------------------------------------------------
  # Check input parameters
  # ---------------------------------------------------------------------
  eps <- 1
  while (1 + eps != 1) eps = eps/2; eps = 2 * eps
  # Profiling option
  options(digits.secs = 3)
  m <- nrow(A)
  n <- ncol(A)
  k <- max(nu, nv)
  interchange <- FALSE
  # Interchange dimensions m, n so that dim(A'A) = min(m, n) when seeking
  # the smallest singular values. This avoids finding zero smallest
  # singular values.
  if (n>m && sigma =="ss") {
    t <- m
    m <- n
    n <- t
    interchange <- TRUE
  }

  # Increase the number of desired signular values by 'adjust' to
  # help convergence. k is re-adjusted as vectors converge--this is
  # only an initial value;
  k_org <- k;
  k <- k + adjust;
  if (k <= 0)  stop ("k must be positive")
  if (k>min(m, n)) stop ("k must be less than min(m, n)+ adjust")
  if (m_b <= 1) stop ("m_b must be greater than 1")
  if (tol<0) stop ("tol must be non-negative")
  if (maxit <= 0) stop ("maxit must be positive")
  if (m_b >= min(n, m)) {
    m_b <- floor(min(n, m)-0.1)
    if (m_b-k-1<0) {
      adjust <- 0
      k <- m_b-1
    }
  }
  if (m_b-k-1<0) m_b <- ceiling(k + 1 + 0.1)
  if (m_b >= min(m, n)) {
    m_b <- floor(min(m, n)-0.1)
    adjust <- 0
    k <- m_b - 1
  }
  if (tol<eps) tol <- eps

  # Allocate memory for W and F:
  W <- matrix(0.0, m, m_b)
  F <- matrix(0.0, n, 1)
  # If starting matrix V is not given then set V to be an
  # (n x 1) matrix of normally distributed random numbers.
  # In any case, allocate V appropriate to problem size:
  if (is.null(V)) {
    V <- matrix(0.0, n, m_b)
    V[, 1] <- rnorm(n)
  }
  else {
    V <- cbind(V, matrix(0.0, n, m_b-ncol(V)))
  }


  # ---------------------------------------------------------------------
  # Initialize local variables
  # ---------------------------------------------------------------------

  B <- NULL                  # Bidiagonal matrix
  Bsz <- NULL                # Size of B
  eps23 <- eps^(2/3)         # Used for Smax/avoids using zero
  I <- NULL                  # Indexing
  J <- NULL                  # Indexing
  iter <- 1                  # Man loop iteration count
  mprod <- 0                 # Number of matrix-vector products
  R_F <- NULL                # 2-norm of residual vector F
  sqrteps <- sqrt(eps)       #
  Smax <- 1                  # Max value of all computed singular values of
  # B est. ||A||_2
  Smin <- NULL               # Min value of all computed singular values of
  # B est. cond(A)
  SVTol <- min(sqrteps, tol)  # Tolerance for singular vector convergence
  S_B <- NULL                # Singular values of B
  U_B <- NULL                # Left singular vectors of B
  V_B <- NULL                # Right singular vectors of B
  V_B_last <- NULL           # last row of modified V_B
  S_B2 <- NULL               # S.V. of [B ||F||]
  U_B2 <- NULL               #
  V_B2 <- NULL               #

  # ---------------------------------------------------------------------
  # Basic functions
  # ---------------------------------------------------------------------

  # Euclidean norm
  norm <- function (x) return(as.numeric(sqrt(crossprod(x))))

  # Orthogonalize vectors Y against vectors X. Y and X must be R matrix
  # objects (they must have a dim attribute).
  # Note: this function unnecessarily copies the contents of Y
  orthog <- function (Y, X) {
    if (dim(X)[2] < dim(Y)[2]) dotY <- crossprod (X, Y)
    else dotY <- t (crossprod(Y, X))
    return (Y - X %*% dotY)
  }

  # Convergence tests
  # Input parameters
  # Bsz            Number of rows of the bidiagonal matrix B
  # tol
  # k_org
  # U_B            Left singular vectors of small matrix B
  # S_B            Singular values of B
  # V_B            Right singular vectors of B
  # residuals
  # k
  # SVTol
  # Smax
  #
  # Output parameter list
  # converged      TRUE/FALSE
  # U_B            Left singular vectors of small matrix B
  # S_B            Singular values of B
  # V_B            Right singular vectors of B
  # k              Number of singular vectors returned
  convtests <- function (Bsz, tol, k_org, U_B, S_B, V_B,
    residuals, k, SVTol, Smax) {
    Len_res <- sum(residuals[1:k_org] < tol * Smax)
    if (Len_res == k_org) {
      return (list(converged = TRUE, U_B = U_B[, 1:k_org, drop = FALSE],
        S_B = S_B[1:k_org, drop = FALSE], V_B = V_B[, 1:k_org, drop = FALSE],
        k = k) )
    }
    #   Not converged yet...
    #   Adjust k to include more vectors as the number of vectors converge.
    Len_res <- sum(residuals[1:k_org] < SVTol * Smax)
    k <- max(k, k_org + Len_res)
    if (k > Bsz -3) k <- Bsz -3
    return (list(converged = FALSE, U_B = U_B, S_B = S_B, V_B = V_B, k = k) )
  }

  # ---------------------------------------------------------------------
  # Main iteration
  # ---------------------------------------------------------------------

  while (iter <= maxit) {

    # ---------------------------------------------------------------------
    # Lanczos bidiagonalization iteration
    # Compute the Lanczos bidiagonal decomposition:
    # AV  = WB
    # A'W = VB + FE'
    # with full reorthogonalization.
    # This routine updates W, V, F, B, mprod
    # ---------------------------------------------------------------------
    j <- 1
    #   Normalize starting vector:
    if (iter == 1) V[, 1] <- V[, 1, drop = FALSE]/norm(V[, 1, drop = FALSE])
    else j <- k + 1

    #   Compute W = AV (the as.matrix fn converts Matrix class objects)
    if (interchange)  W[, j] <- t (as.matrix(crossprod (V[, j, drop = FALSE], A)))
    else              W[, j] <- as.matrix(A %*% V[, j, drop = FALSE])
    mprod <- mprod + 1

    #   Orthogonalize
    if (iter != 1) {
      W[, j] <- orthog (W[, j, drop = FALSE], W[, 1:j-1, drop = FALSE])
    }

    S <- norm(W[, j, drop = FALSE])
    #   Check for linearly-dependent vectors
    if ((S < SVTol) && (j == 1)) stop ("Starting vector near the null space")
    if (S < SVTol) {
      W[, j] <- rnorm(nrow(W))
      W[, j] <- orthog(W[, j, drop = FALSE], W[, 1:j-1, drop = FALSE])
      W[, j] <- W[, j, drop = FALSE]/norm(W[, j, drop = FALSE])
      S <- 0
    }
    else W[, j] <- W[, j, drop = FALSE]/S

    #   Lanczos process
    while (j <= m_b) {
      if (interchange) F <- as.matrix(A %*% W[, j, drop = FALSE])
      else F <- t(as.matrix(crossprod(W[, j, drop = FALSE], A)))

      mprod <- mprod + 1
      F <- F - S * V[, j, drop = FALSE]
      #     Orthogonalize
      F <- orthog(F, V[, 1:j, drop = FALSE])

      if (j + 1 <= m_b) {
        R <- norm(F)
        #       Check for linear dependence
        if (R <= SVTol) {
          F <- matrix(rnorm(dim(V)[1]), dim(V)[1], 1)
          F <- orthog(F, V[, 1:j, drop = FALSE])
          V[, j + 1] <- F/norm(F)
          R <- 0
        }
        else V[, j + 1] <- F/R

        #       Compute block diagonal matrix
        if (is.null(B)) B <- cbind(S, R)
        else            B <- rbind(cbind(B, 0), c(rep(0, j-1), S, R))

        if (interchange) W[, j + 1] <- t (as.matrix(crossprod (V[, j + 1,
          drop = FALSE], A)))
        else             W[, j + 1] <- as.matrix(A %*% V[, j + 1, drop = FALSE])
        mprod <- mprod + 1

        #       One step of the classical Gram-Schmidt process
        W[, j + 1] <- W[, j + 1, drop = FALSE] - W[, j, drop = FALSE]* R

        #       Full reorthogonalization of W
        if (iter == 1 || reorth == 2)
          W[, j + 1] <- orthog(W[, j + 1, drop = FALSE], W[, 1:j, drop = FALSE])
        S <- norm(W[, j + 1, drop = FALSE])
        if (S <= SVTol) {
          W[, j + 1] <- rnorm(nrow(W))
          W[, j + 1] <- orthog(W[, j + 1, drop = FALSE], W[, 1:j, drop = FALSE])
          W[, j + 1] <- W[, j + 1, drop = FALSE]/norm(W[, j + 1, drop = FALSE])
          S <- 0
        }
        else W[, j + 1] <- W[, j + 1, drop = FALSE]/S
      }
      else {
        #       Add a last block to matrix B
        B <- rbind(B, c(rep(0, j - 1), S))
      }
      j <- j + 1
    }
    #cat ("iter = ", iter," j = ", j-1, "mprod = ", mprod,"\n", file = stderr())
    # ---------------------------------------------------------------------
    # (End of the Lanczos bidiagonalization part)
    # ---------------------------------------------------------------------

    Bsz <- nrow(B)
    R_F <- norm(F)
    F <- F/R_F
    #   Compute singular triplets of B. Expect svd to return s.v.s in order
    #   from largest to smallest.
    Bsvd <- svd(B)

    #   Estimate ||A|| using the largest singular value over all iterations
    #   and estimate the cond(A) using approximations to the largest and
    #   smallest singular values. If a small singular value is less than sqrteps
    #   use only Ritz vectors to augment and require two-sided reorthogonalization.
    if (iter == 1) {
      Smax <- Bsvd$d[1]
      Smin <- Bsvd$d[Bsz]
    }
    else {
      Smax <- max(Smax, Bsvd$d[1])
      Smin <- min(Smin, Bsvd$d[Bsz])
    }
    Smax <- max(eps23, Smax)
    if ((Smin/Smax < sqrteps) && reorth <2) {
      warning ("The matrix is ill-conditioned.",
        "Each basis will be reorthogonalized.")
      reorth <- 2
      aug <- "ritz"
    }

    #   Re-order the singular values accordingly.
    if (sigma == "ss") {
      jj <- seq (ncol (Bsvd$u), 1, by = -1)
      Bsvd$u <- Bsvd$u[, jj]
      Bsvd$d <- Bsvd$d[jj]
      Bsvd$v <- Bsvd$v[, jj]
    }

    #   Compute the residuals
    R <- R_F %*% Bsvd$u[Bsz,, drop = FALSE]
    #   Check for convergence
    ct <- convtests(Bsz, tol, k_org, Bsvd$u, Bsvd$d, Bsvd$v, abs(R), k, SVTol,
      Smax)

    #   If all desired singular values converged, then exit main loop
    if (ct$converged) break

    if (iter >= maxit) break

    #   Compute starting vectors & 1st block of B[1:k, 1:(k + 1), drop = FALSE]
    if (aug == "harm") {
      #     Update the SVD of B to be the svd of [B ||F||E_m]
      Bsvd2 <- svd (cbind (diag (Bsvd$d), t (R)))
      if (sigma == "ss") {
        jj <- seq (ncol (Bsvd2$u), 1, by =-1)
        Bsvd2$u <- Bsvd2$u[, jj]
        Bsvd2$d <- Bsvd2$d[, jj]
        Bsvd2$v <- Bsvd2$v[, jj]
      }
      Bsvd$d <- Bsvd2$d
      Bsvd$u <- Bsvd$u %*% Bsvd2$u
      Bsvd$v <- cbind (rbind (Bsvd$v, rep (0, Bsz)),
        c (rep (0, Bsz), 1)) %*% Bsvd2$v
      V_B_last <- Bsvd$v [Bsz + 1, , drop = FALSE]
      s <- R_F %*% solve (B, cbind (c (rep(0, Bsz - 1), 1)))
      Bsvd$v <- Bsvd$v[1:Bsz, , drop = FALSE] + s * Bsvd$v[Bsz + 1, ,
        drop = FALSE]

      qrv <- qr (cbind ( rbind (Bsvd$v[, 1:k], 0), rbind (-s, 1)))
      Bsvd$v[, 1:k + 1] <- cbind (Bsvd$v, F) %*% qr.Q (qrv)

      #     Update and compute the k x k + 1 part of B
      UT <- R[1:k + 1, 1:k, drop = FALSE] + R[, k + 1, drop = FALSE] %*%
        t(V_B_last)
      B <- diag (Bsvd$d[1:k]) %*% (UT * upper.tri (UT, diag = TRUE))
    }
    else {
      #     Use the Ritz vectors
      V[, 1:(k + dim(F)[2])] <- cbind(V[, 1:(dim(Bsvd$v)[1]), drop = FALSE] %*%
          Bsvd$v[, 1:k, drop = FALSE], F)
      B <- cbind( diag(Bsvd$d[1:k, drop = FALSE]), R[1:k, drop = FALSE])
    }

    #   Update the left approximate singular vectors
    W[, 1:k] <- W[, 1:(dim(Bsvd$u)[1]), drop = FALSE] %*% Bsvd$u[, 1:k,
      drop = FALSE]

    iter <- iter + 1
  }
  # ---------------------------------------------------------------------
  # End of the main iteration loop
  # ---------------------------------------------------------------------

  # ---------------------------------------------------------------------
  # Output results
  # ---------------------------------------------------------------------

  #DPT. OF DIRTY TRICKS: modified names to conform with those returned by eigen
  #to save some keystrokes/assignments in orthoDesign() ...
  return ( list(values = Bsvd$d[1:k_org], vectors = W[, 1:(dim(Bsvd$u)[1]),
    drop = FALSE] %*% Bsvd$u[, 1:k_org, drop = FALSE],
    v = V[, 1:(dim(Bsvd$v)[1]), drop = FALSE] %*% Bsvd$v[, 1:k_org, drop = FALSE],
    iter = iter))
}

# like cat, but repects width
wrapCat <- function(...) {
  str <- as.vector(list(...))
  fmt <- strwrap(str, width = getOption("width"), exdent = 2)
  fmt <- paste(fmt, collapse ="\n")
  cat(fmt)
}
