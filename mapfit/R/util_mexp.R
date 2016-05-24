#########

mexp.options <- function() {
    list(
        poi.eps = sqrt(.Machine$double.eps),
        unif.fact = 1.01,
        unif.tol = sqrt(.Machine$double.eps),
        pade.eps = sqrt(.Machine$double.eps),
        arnoldi.iter = 1,
        arnoldi.tol = sqrt(.Machine$double.eps)
        )
}

mexp.pade <- function(transpose, A, eps = sqrt(.Machine$double.eps)) {
  if (is.matrix(A)) {
    A <- as(A, "dgeMatrix")
  }
  if (transpose == TRUE) {
    .Call(mexp_pade, dim(A)[1L], Matrix::t(A), eps)
  } else {
    .Call(mexp_pade, dim(A)[1L], A, eps)
  }
}

mexp.uniform <- function(transpose, A, t = 1.0, eps = sqrt(.Machine$double.eps),
  ufact = 1.01, atol = sqrt(.Machine$double.eps)) {
  if (is.matrix(A)) {
    A <- as(A, "dgeMatrix")
  }
  if (transpose == TRUE) {
    .Call(mexp_unif, dim(A)[1L], t(A), t, eps, ufact, atol)
  } else {
    .Call(mexp_unif, dim(A)[1L], A, t, eps, ufact, atol)
  }
}

mexp.unifvec <- function(transpose, A, x, t = 1.0, eps = sqrt(.Machine$double.eps), 
  ufact = 1.01, atol = sqrt(.Machine$double.eps), class = "CsparseMatrix") {
  if (is.matrix(A)) {
    A <- as(A, class)
  }
  .Call(mexp_unifvec, transpose, dim(A)[1L], A, x, t, eps, ufact, atol)
}

mexp.unifseq <- function(transpose, A, x, t, eps = sqrt(.Machine$double.eps), 
  ufact = 1.01, atol = sqrt(.Machine$double.eps), class = "CsparseMatrix") {
  if (is.matrix(A)) {
    A <- as(A, class)
  }
  result <- .Call(mexp_unifseq, transpose, dim(A)[1L], A, x, diff(t), eps, ufact, atol)
  list(time=t, result=result)
}

mexp.kryvec <- function(transpose, A, x, t = 1.0, ksub = 10,
  ite = 1, tol = sqrt(.Machine$double.eps), eps = sqrt(.Machine$double.eps), class = "CsparseMatrix") {
  if (is.matrix(A)) {
    A <- as(A, class)
  }
  res <- .Call(mexp_kryvec, transpose, dim(A)[1L], A, x, t, ksub, ite, tol, eps)
  list(value=res[[1L]], err=res[[2L]])
}

mexp <- function(A, x, t = 1.0, transpose = FALSE,
  method = c("automatic", "pade", "uniformization", "krylov"), ksub = 30,
  options = list(), ...) {
    method <- match.arg(method)
    call <- match.call()

    con <- mexp.options()
    nmsC <- names(con)
    con[(namc <- names(options))] <- options
    if (length(noNms <- namc[!namc %in% nmsC]))
        warning("unknown names in option: ", paste(noNms, collapse = ", "))

    if (missing(x)) {
      if (method == "automatic") {
        method <- "pade"
      }
      switch(method,
        "pade" = mexp.pade(transpose, as.matrix(A*t), con$pade.eps, ...),
        "uniformization" = mexp.uniform(transpose, as.matrix(A), t=t, eps=con$poi.eps,
          ufact=con$unif.fact, atol=con$unif.tol, ...),
        "krylov" = stop("Krylov method cannot be applied to computing the full matrix exponential."),
        stop("unknown method: ", method)
        )
    } else if (length(t) == 1) {
      if (method == "automatic") {
        if (is.matrix(A) || is(A, "dgeMatrix")) {
          method <- "pade"
        } else if (is(A, "dgCMatrix") || is(A, "dgRMatrix") || is(A, "dgTMatrix")) {
          method <- "uniformization"
        } else {
          method <- "uniformization"
          A <- as(A, "dgTMatrix")
        }
      }
      if (method == "pade") {
        if (!is(A,"dgeMatrix")) {
          A <- as(A, "dgeMatrix")
        }
      }
      switch(method,
        "pade" = c(mexp.pade(transpose, A*t, eps=con$pade.eps) %*% x, ...),
        "uniformization" = mexp.unifvec(transpose, A, x, t=t,
          eps=con$poi.eps, ufact=con$unif.fact, atol=con$unif.tol, ...),
        "krylov" = mexp.kryvec(transpose, A, x, t=t, ksub=ksub,
          eps=con$pade.eps, ite=con$arnoldi.iter, tol=con$arnoldi.tol, ...),
        stop("unknown method: ", method)
      )
    } else {
      method <- "uniformization"
      switch(method,
        "uniformization" = mexp.unifseq(transpose, A, x, t,
          eps=con$poi.eps, ufact=con$unif.fact, atol=con$unif.tol),
        stop("unknown method: ", method)
      )
    }
}

mpow <- function(A, m, transpose = FALSE) {
  if (is.matrix(A)) {
    A <- as(A, "dgeMatrix")
  }
  if (transpose == TRUE) {
    .Call(mpow_mat, t(A), m)
  } else {
    .Call(mpow_mat, A, m)
  }
}

msolve <- function(alpha, A, x, transpose = FALSE, maxiter = 10000,
  eps = sqrt(.Machine$double.eps)) {
  if (is.matrix(A)) {
    A <- as(A, "dgeMatrix")
  }
  .Call(marsolve, transpose, alpha, A, x, maxiter, eps)
}

ctmc.st <- function(Q, maxiter = 2000, eps = sqrt(.Machine$double.eps)) {
  if (is.matrix(Q)) {
    Q <- as(Q, "dgeMatrix")
  }
  .Call(ctmc_st, Q, maxiter, eps)
}

