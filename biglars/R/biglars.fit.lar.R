"biglars.fit.lar" <- 
function(R, Qty, removeColumns = TRUE, eps = sqrt(.Machine$double.eps), ...)
{
  vecnorm <- function(x)
    sqrt(sum(x * x))

  formXtr <- function(R, Qty, b, step = 0, d = NULL)
  {
    if(step && !is.null(d))
      b <- b + step * d
    drop(crossprod(R, (Qty - R %*% b)))
  }
  Rdrop <- function(R, eps)
  {
    rmag <- abs(diag(R))
    rmag < eps * (1 + max(rmag))
  }
  dimR <- nrow(R)
  if(ncol(R) != dimR)
    stop("R should be square")
  if(any(as.logical(R[row(R) > col(R)])))
    stop("R should be upper triangular")
  RtR <- crossprod(R)
  ADD <- DROP <- rep(0, dimR)
  # always include the intercept
  colNorms <- apply(R, 2, vecnorm)
  if(removeColumns) {
    BAD <- colNorms < eps * (1 + max(colNorms))
    if(any(BAD))
      DROP[BAD] <- -1
  }
  b <- matrix(0, dimR, dimR)
  storage.mode(R) <- "double"
  k <- 1
  ACTIVE <- NULL
  while(TRUE) {
    gvec <- formXtr(R, Qty, b[k,  ])
    cvec <- gvec/colNorms
    cmax <- max(abs(cvec[!DROP]))
    NEW <- abs(abs(cvec) - cmax) < eps * (1 + cmax) & !DROP
    NEW <- setdiff(which(NEW), ACTIVE)
    ACTIVE <- sort(c(ACTIVE, NEW))
    ADD[NEW] <- k
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
      PACKAGE = "biglars")[c('R', 'y')]
    restoreRy$R <- restoreRy$R[1:nACTIVE,  , drop = FALSE]
    restoreRy$R[row(restoreRy$R) > col(restoreRy$R)] <- 0
    restoreRy$y <- restoreRy$y[1:nACTIVE]
    same <- FALSE
    if(any(rDrop <- Rdrop(restoreRy$R, eps))) {
      DROP[ACTIVE[rDrop]] <-  - k
      BAD <- which(DROP ==  - k)
      if(length(BAD) != length(NEW) || !all(NEW == BAD)) {
        ACTIVE <- which(ADD & !DROP)
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
          PACKAGE = "biglars")[c('R', 'y')]
        restoreRy$R <- restoreRy$R[1:nACTIVE,  , drop = FALSE]
        restoreRy$R[row(restoreRy$R) > col(restoreRy$R)] <- 0
        restoreRy$y <- restoreRy$y[1:nACTIVE]
      }
      else same <- TRUE
    }
    if(!same) {
      d <- rep(0, dimR)
      ##
      ## EQUIVALENT
      ##
      ##   z <- solve(restoreRy$R[1:nACTIVE, 1:nACTIVE, drop = FALSE], 
      ##              restoreRy$y[1:nACTIVE] - restoreRy$R %*% b[k, ACTIVE])
      ##      d[ACTIVE] <- z
      ##
      z <- solve(restoreRy$R[1:nACTIVE, 1:nACTIVE, drop = FALSE], restoreRy$y[1:
        nACTIVE])
      d[ACTIVE] <- z - b[k, ACTIVE]
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
      if(FALSE) {
        gvec <-  - crossprod(R, R %*% d)
        dvec <- gvec/colNorms
        A <- 0
        B <- 1
        gamma <- 1/2
        FREE <- which(!ADD & !DROP)
        ACT <- ACTIVE
        while(TRUE) {
          v <- gvec + gamma * dvec
          print(c(A, gamma, B))
          print(cbind(as.numeric(ADD | DROP), v, gvec, dvec))
          vmag <- max(abs(v))
          vFREE <- max(abs(v[FREE]))
          vACT <- max(abs(v[ACT]))
          print(c(vACT, vFREE, vmag))
          if(abs(vFREE - vACT) < eps * vmag)
            break
          if(vFREE < vACT) {
            A <- A + (B - A)/2
          }
          else {
            B <- B - (B - A)/2
          }
          gamma <- (A + B)/2
          if(abs(B - A) < eps)
            stop("STOP")
        }
        print(gamma)
      }
      sact <- mean((svec[ACTIVE]))
      dact <- mean((del[ACTIVE]))
      EPS <- eps * (1 + max(svec))
      bmax <- 1 + max(abs(b[k,  ]))
      dmax <- 1 + max(abs(d))
      gamhat <- 1
      FREE <- which(!ADD & !DROP)
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
      b[k, ACTIVE] <- b[k, ACTIVE] + gamhat * d[ACTIVE]
      if((nACTIVE + sum(as.logical(DROP))) >= dimR)
        break
      b[k + 1,  ] <- b[k,  ]
      k <- k + 1
    }
  }
  moves <- cbind(Var = 1:dimR, Stage = ADD)
  if(any(as.logical(DROP))) {
    warning("columns dropped due to ill-conditioning: ", paste(which(DROP != 0),
      collapse = " "))
    moves <- rbind(moves, cbind(Var =  - (1:dimR), Stage = abs(DROP)))
  }
  moves <- moves[moves[, "Stage"] > 0,  , drop = FALSE]
  moves <- moves[order(moves[, "Stage"]),  , drop = FALSE]
  list(coef = b[1:k,  ], moves = moves, keep = which(as.logical(ADD)))
}

