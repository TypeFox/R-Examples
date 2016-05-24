"biglars.fit.stepwise" <- 
function(R, Qty, removeColumns = TRUE, eps = sqrt(.Machine$double.eps), ...)
{
  if(is.R()) {
    vecnorm <- function(x)
    sqrt(sum(x * x))
  }
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
  ADD[1] <- 1
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
    ADD[NEW] <- k
    ACTIVE <- sort(c(ACTIVE, NEW))
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
          PACKAGE = "biglars")[c("R", "y")]
        restoreRy$R <- restoreRy$R[1:nACTIVE,  , drop = FALSE]
        restoreRy$R[row(restoreRy$R) > col(restoreRy$R)] <- 0
        restoreRy$y <- restoreRy$y[1:nACTIVE]
      }
      else same <- TRUE
    }
    if(!same) {
      ##      z <- solve(restoreRy$R[1:nACTIVE, 1:nACTIVE, drop = FALSE], 
      ##                 restoreRy$y[1:nACTIVE] - restoreRy$R %*% b[k, ACTIVE])
      ##     d[ACTIVE] <- z
      z <- solve(restoreRy$R[1:nACTIVE, 1:nACTIVE, drop = FALSE], restoreRy$y[1:
        nACTIVE])
      b[k, ACTIVE] <- z
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

