"biglars.fit.lasso" <- 
function(R, Qty, removeColumns = TRUE, eps = sqrt(.Machine$double.eps), 
  maxStages = NULL)
{
  vecnorm <- function(x)
    sqrt(sum(x * x))

  formXtr <- function(R, Qty, b, step = 0, d = NULL, ACTIVE = NULL)
  {
    if(step && !is.null(d)) {
      crossprod(R, (Qty - R %*% (b + step * d)))
    }
    else {
      crossprod(R, (Qty - R %*% b))
    }
  }
  Rdrop <- function(R, eps)
  {
    rmag <- abs(diag(R))
    rmag < eps * (1 + max(rmag))
  }
  updateDROP <- function(DROP, MOVES)
  {
    ACTIVE <- apply(MOVES, 2, function(x)
    sum(sign(x) > 0))
    for(k in unique(DROP)) {
      if(k == 1)
        next
      K <- DROP == k
      ACTIVEk <- apply(MOVES[1:k,  , drop = FALSE], 2, function(x)
      sum(sign(x)) > 0)
      if(any(ACTIVEk & !ACTIVE))
        DROP[K] <- 0
    }
    DROP
  }
  changeSign <- function(b, d, eps)
  {
    beps <- eps * (1 + max(abs(b)))
    deps <- eps * (1 + max(abs(d)))
    b[abs(b) < beps] <- 0
    d[abs(d) < deps] <- 0
    S <- rep(Inf, length(b))
    if(any(I <- (b & sign(b) != sign(d))))
      S[I] <-  - b[I]/d[I]
    S
  }
  lassoCond <- function(b, g, ACTIVE, eps)
  {
    beps <- eps * (1 + max(abs(b)))
    gmax <- max(abs(g))
    geps <- eps * (1 + gmax)
    b[abs(b) < beps] <- g[abs(g) < geps] <- 0
    ACT <- b & sign(b) == sign(g)
    # active set
    BIG <- abs(abs(g) - gmax) < geps
    all(BIG[ACTIVE])
  }
  dimR <- nrow(R)
  if(ncol(R) != dimR)
    stop("R should be square")
  if(any(as.logical(R[row(R) > col(R)])))
    stop("R should be upper triangular")
  RtR <- crossprod(R)
  # always include the intercept
  colNorms <- apply(R, 2, vecnorm)
  if(removeColumns) {
    REMOVE <- colNorms < eps * (1 + max(colNorms))
  }
  else REMOVE <- rep(FALSE, dimR)
  DROP <- rep(FALSE, dimR)
  DROP[REMOVE] <- 1
  if(is.null(maxStages))
    maxStages <- 2 * dimR
  CURRENT <- rep(0, dimR)
  b <- MOVES <- matrix(0., maxStages, dimR)
  beta <- rep(0, dimR)
  # coef for intercept only
  # Find variables with the highest correlation with y
  # k = step number; begin step 2 (step 1 is preliminary)
  k <- 1
  storage.mode(R) <- "double"
  gvec <- formXtr(R, Qty, b[k,  ])
  cvec <- gvec/colNorms
  cmax <- max(abs(cvec[!REMOVE]))
  NEW <- abs(abs(cvec) - cmax) < eps * (1 + cmax) & !REMOVE
  ACTIVE <- which(NEW)
  MOVES[k, NEW] <- 1
  while(TRUE) {
    nACTIVE <- length(ACTIVE)
    restoreRy <- .Fortran("rstrup",
      as.integer(dimR),
      as.integer(nACTIVE),
      R = R[, ACTIVE, drop = FALSE],
      as.integer(dimR),
      as.integer(ACTIVE),
      double(dimR),
      double(dimR),
      y = as.double(Qty),
      PACKAGE = "biglars")[c("R", "y")]
    restoreRy$R <- restoreRy$R[1:nACTIVE,  , drop = FALSE]
    restoreRy$R[row(restoreRy$R) > col(restoreRy$R)] <- 0
    restoreRy$y <- restoreRy$y[1:nACTIVE]
    same <- FALSE
    BAD <- NULL
    if(any(rDrop <- Rdrop(restoreRy$R, eps))) {
      MOVES[k, ACTIVE[rDrop]] <-  - Inf
      BAD <- which(MOVES[k,  ] ==  - Inf)
      DROP[BAD] <- k
      if(length(BAD) != length(NEW) || !all(NEW == BAD)) {
        ACTIVE <- setdiff(ACTIVE, BAD)
        nACTIVE <- length(ACTIVE)
        restoreRy <- .Fortran("rstrup",
          as.integer(dimR),
          as.integer(nACTIVE),
          R = R[, ACTIVE, drop = FALSE],
          as.integer(dimR),
          as.integer(ACTIVE),
          double(dimR),
          double(dimR),
          y = as.double(Qty),
          PACKAGE = "biglars")[c("R", "y")]
        restoreRy$R <- restoreRy$R[1:nACTIVE,  , drop = FALSE]
        restoreRy$R[row(restoreRy$R) > col(restoreRy$R)] <- 0
        restoreRy$y <- restoreRy$y[1:nACTIVE]
      }
      else same <- TRUE
    }
    if(!same) {
      z <- solve(restoreRy$R[1:nACTIVE, 1:nACTIVE, drop = FALSE], restoreRy$y[1:
        nACTIVE])
      d <- rep(0, dimR)
      d[ACTIVE] <- z - beta[ACTIVE]
      # OLS regression coefficients, for res against active x's
      # Move in direction of OLS solution, but only as long as
      # correlation with partial residuals is largest for active
      # variables.  Compute gamhat = the fraction of the
      # distance to the full LS solution to go.  At gamhat = 1, 
      # the correlation of active variables with the partial
      # residuals would be zero.
      s <- sign(gvec)
      s[!s] <- 1
      del <- (s * crossprod(R, R %*% d))/colNorms
      svec <- abs(cvec)
      sact <- mean(svec[ACTIVE])
      dact <- mean(del[ACTIVE])
      EPS <- eps * (1 + max(svec))
      gamhat <- 1
      FREE <- setdiff(which(!REMOVE), c(BAD, ACTIVE))
      bmax <- 1 + max(abs(beta))
      dmax <- 1 + max(abs(d))
      for(i in FREE) {
        top <- sact - svec[i]
        bot <- dact - del[i]
        skip <- abs(top) > abs(bot) | abs(bot) < EPS | sign(top) != sign(bot)
        skip <- skip || abs(top) * dmax < eps * abs(bot) * bmax
        if(!skip)
          gamhat <- min(gamhat, top/bot)
        top <- svec[i] + sact
        bot <- del[i] + dact
        skip <- abs(top) > abs(bot) | abs(bot) < EPS | sign(top) != sign(bot)
        skip <- skip || abs(top) * dmax < eps * abs(bot) * bmax
        if(!skip)
          gamhat <- min(gamhat, top/bot)
      }
      step <- min(changeSign(beta[ACTIVE], d[ACTIVE], eps))
      if(step > gamhat) {
        beta <- beta + gamhat * d
        ## the drops have to be handled properly using MOVES
        gvec <- formXtr(R, Qty, beta)
        cvec <- gvec/colNorms
        if(lassoCond(beta, cvec, ACTIVE, eps)) {
          b[k,  ] <- beta
          if(nACTIVE + sum(as.logical(DROP)) >= dimR)
            break
          k <- k + 1
        }
        else stop("STOP")
        cmax <- max(abs(cvec[!REMOVE]))
        OLD <- ACTIVE
        NEW <- abs(abs(cvec) - cmax) < eps * (1 + cmax) & !REMOVE
        NEW <- setdiff(which(NEW), ACTIVE)
        MOVES[k, NEW] <- 1
        ACTIVE <- sort(c(ACTIVE, NEW))
        if(is.null(setdiff(OLD, ACTIVE)) && is.null(setdiff(ACTIVE, OLD)))
          break
      }
      else {
        beta <- beta + step * d
        Z <- abs(beta) < eps * (1 + max(abs(beta)))
        Z <- intersect(ACTIVE, which(Z))
        MOVES[k, Z] <- -1
        gvec <- formXtr(R, Qty, beta)
        cvec <- gvec/colNorms
        if(lassoCond(beta, cvec, ACTIVE, eps)) {
          b[k,  ] <- beta
          k <- k + 1
        }
        else stop("STOP")
        cmax <- max(abs(cvec[!REMOVE]))
        OLD <- ACTIVE
        NEW <- abs(abs(cvec) - cmax) < eps * (1 + cmax) & !REMOVE
        NEW <- setdiff(which(NEW), ACTIVE)
        ## the drops have to be handled properly using MOVES
        ACTIVE <- sort(c(ACTIVE, NEW))
        ACTIVE <- setdiff(ACTIVE, Z)
        if(is.null(setdiff(OLD, ACTIVE)) && is.null(setdiff(ACTIVE, OLD)))
          break
        MOVES[k, NEW] <- k
      }
      DROP <- updateDROP(DROP, MOVES)
    }
  }
  MOVES <- sign(MOVES[1:k,  , drop = FALSE])
  STAGE <- as.vector(t(sweep(abs(MOVES), MARGIN = 1, FUN = "*", STATS = 1:k)))
  MOVES <- as.vector(t(sweep(MOVES, MARGIN = 2, FUN = "*", STATS = 1:dimR)))
  l <- MOVES != 0
  moves <- cbind(Var = MOVES[l], Stage = STAGE[l])
  list(coef = b[1:k,  ], moves = moves)
}

