# s_sd <- function(x, mu.too = FALSE, ...) {
#   c(if (mu.too) mean(x), sd(x, ...))
# }

robcor <- function(x, y = NULL, method = c("ssd", "quadrant", "mcd"), partial = FALSE, post = "psdcor", scaler = "s_FastQn", regress = "lmrob") {
  method <- match.arg(method)
  scaler <- match.fun(scaler)
  if (!is.null(post))
    post <- match.fun(post)
  if (is.data.frame(y)) 
    y <- as.matrix(y)
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!is.matrix(x) && is.null(y)) 
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  if (!is.null(y)) {
    if (!(is.numeric(y) || is.logical(y)))
      stop("'y' must be numeric")
    stopifnot(is.atomic(y))
    if (method == "mcd")
      stop(paste("'y' must be NULL for", method))
    if (partial)
      stop("'y' must be NULL for partial correlations")
  }

  if (partial && any(method == c("mcd")))
    stop(paste("partial correlations unavailable for", method))

  if (method == "mcd") {
    ret <- covMcd(x, cor = TRUE)$cor
  } else {
    ret <- .PairwiseCorrelation(x, y, method, partial, post, scaler, regress)
  }
  ret
}

.PairwiseCorrelation <- function(x, y, method, partial, post, scaler, regress) {
  if (method == "ssd") {
    .Standardize <- function(x) {
      musigma <- apply(x, 2, function(x) scaler(na.omit(x), mu.too=TRUE))
      scale(x, musigma[1,], musigma[2,])
    }
    .PairCorrelation <- function(nx, ny) {
      s2u <- scaler(nx + ny)^2
      s2v <- scaler(nx - ny)^2
      (s2u - s2v) / (s2u + s2v)
    }
  } else if (method == "quadrant") {
    .Standardize <- function(x) {
      center <- apply(x, 2, median, na.rm = TRUE)
      sign(sweep(x, 2, center))
    }
    .PairCorrelation <- function(nx, ny) {
      n <- nx * ny
      sin(sum(n) / sum(abs(n)) * pi / 2)
    }
  } else {
    stop("unimplemented method")
  }

  if (is.null(y)) {
    ncy <- ncx <- ncol(x)
    if (ncx == 0) 
      stop("'x' is empty")
    if (!is.null(.Standardize))
      x <- .Standardize(x)

    if (partial) {
      if (regress == "lmrob") {
        require(robustbase)
      }
      regress <- match.fun(regress)
      r <- .PartialToCorrelation(.BuildPartial(x, .Standardize, .PairCorrelation, regress))
    } else {
      r <- diag(ncx)
      for (i in seq_len(ncx - 1)) {
        for (j in i + seq_len(ncx - i)) {
          x2 <- x[, i]
          y2 <- x[, j]
          ok <- complete.cases(x2, y2)
          x2 <- x2[ok]
          y2 <- y2[ok]
          r[i, j] <- r[j, i] <- if (any(ok)) {
            .PairCorrelation(x2, y2)
          } else {
            NA
          }
        }
      }
    }
    if (is.function(post))
      r <- post(r)
    rownames(r) <- colnames(x)
    colnames(r) <- colnames(x)
  } else {
    if (length(x) == 0L || length(y) == 0L) 
        stop("both 'x' and 'y' must be non-empty")
    matrixResult <- is.matrix(x) || is.matrix(y)
    if (!is.matrix(x)) 
        x <- matrix(x, ncol = 1L)
    if (!is.matrix(y)) 
        y <- matrix(y, ncol = 1L)
    if (!is.null(.Standardize)) {
      x <- .Standardize(x)
      y <- .Standardize(y)
    }
    ncx <- ncol(x)
    ncy <- ncol(y)
    r <- matrix(0, nrow = ncx, ncol = ncy)
    for (i in seq_len(ncx)) {
      for (j in seq_len(ncy)) {
        x2 <- x[, i]
        y2 <- y[, j]
        ok <- complete.cases(x2, y2)
        x2 <- x2[ok]
        y2 <- y2[ok]
        r[i, j] <- if (any(ok)) {
          .PairCorrelation(x2, y2)
        } else {
          NA
        }
      }
    }
    if (is.function(post))
      r <- post(r)
    rownames(r) <- colnames(x)
    colnames(r) <- colnames(y)
    if (!matrixResult)
      r <- drop(r)
  }
  r
}

.BuildPartial <- function(data, fStandardize, fCorrelate, fRegression) {
  stopifnot(is.matrix(data))

  p <- ncol(data)
  result <- diag(p)

  data <- as.data.frame(data)
  colnames(data) <- paste0("x", seq_len(p))

  for (i in seq_len(p - 1)) {
    if (i == 1) {
      for (j in i + seq_len(p - i)) {
        result[i, j] <- result[j, i] <- fCorrelate(data[,i], data[,j])
      }      
    } else {
      vars <- paste0("x", seq_len(i - 1))
      f.i <- as.formula(paste(paste0("x", i), "~", vars))
      res.i <- residuals(fRegression(f.i, data))
      res.i <- as.vector(fStandardize(as.matrix(res.i)))
      for (j in i + seq_len(p - i)) {
        f.j <- as.formula(paste(paste0("x", j), "~", vars))
        res.j <- residuals(fRegression(f.j, data))
        res.j <- as.vector(fStandardize(as.matrix(res.j)))
        result[i, j] <- result[j, i] <- fCorrelate(res.i, res.j)
      }
    }
  }
  return(result)
}

.PartialToCorrelation <- function(mtx) {
  stopifnot(nrow(mtx) == ncol(mtx))
  stopifnot(isSymmetric(mtx))
  p <- ncol(mtx)
  result <- diag(p)
  for (i in seq_len(p - 1)) {
    for (j in i + seq_len(p - i)) {
      r <- mtx[i, j]
      for (n in seq_len(i - 1)) {
        a <- mtx[i, i - n]
        b <- mtx[i - n, j]
        r <- a * b + r * sqrt(1 - a^2) * sqrt(1 - b^2)
      }
      result[i, j] <- result[j, i] <- r
    }
  }
  result
}
