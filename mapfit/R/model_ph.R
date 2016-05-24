## general PH

ph <- function(size, alpha, Q, xi, class="CsparseMatrix") {
  if (missing(size)) {
    if (missing(alpha) || missing(Q) || missing(xi)) {
      stop("alpha, Q and xi are needed.")
    } else {
      size <- length(alpha)
    }
  } else {
    if (!missing(alpha) && !missing(Q) && !missing(xi)) {
      warning("size is ignored.")
      size <- length(alpha)
    } else {
      if (!missing(alpha) || !missing(Q) || !missing(xi)) {
        warning("alpha, Q and xi are ignored.")
      }
      alpha <- rep(1.0/size, size)
      Q <- matrix(1.0, size, size)
      diag(Q) <- rep(-size, size)
      xi <- rep(1.0, size)
    }
  }
  if (!is(Q, class)) {
    Q <- as(Q, class)
  }
  if (missing(xi)) {
    xi = -apply(Q, 1, sum)
  }
  zero <- 1.0e-8
  df <- sum(abs(alpha) > zero) - 1 + sum(abs(Q) > zero) +
    sum(abs(xi) > zero) - size + sum(abs(Matrix::diag(Q)) < zero)
  new("ph", size=size, alpha=alpha, Q=Q, xi=xi, df=df)
}

ph.bidiag <- function(size, class="CsparseMatrix") {
  if (size <= 1) {
    ph(size, class=class)
  } else {
  alpha <- rep(1/size,size)
  xi <- rep(0, size)
  Q <- matrix(0, size, size)
  for (i in 1:(size-1)) {
    Q[i,i] <- -1
    Q[i,i+1] <- 1
  }
  Q[size,size] <- -1
  xi[size] <- 1
  ph(alpha=alpha, Q=Q, xi=xi, class=class)
  }
}

ph.tridiag <- function(size, class="CsparseMatrix") {
  if (size <= 2) {
    ph(size, class=class)
  } else {
  alpha <- rep(1/size,size)
  xi <- rep(0, size)
  Q <- matrix(0, size, size)
  Q[1,1] <- -1
  Q[1,2] <- 1
  for (i in 2:(size-1)) {
    Q[i,i] <- -2
    Q[i,i-1] <- 1
    Q[i,i+1] <- 1
  }
  Q[size,size-1] <- 1
  Q[size,size] <- -2
  xi[size] <- 1
  ph(alpha=alpha, Q=Q, xi=xi, class=class)
  }
}

setMethod("ph.moment", signature(ph = "ph"),
  function(k, ph, ...) {
  tmp <- ph@alpha
  tmp2 <- 1.0
  res <- numeric(0)
  for (i in 1:k) {
    tmp <- msolve(alpha=1.0, A=-as.matrix(ph@Q), x=tmp, transpose=TRUE)
    tmp2 <- tmp2 * i
    res <- c(res, tmp2 * sum(tmp))
  }
  res
}
)

ph.mean <- function(ph) ph.moment(1, ph)

ph.var <- function(ph) {
  res <- ph.moment(2, ph)
  res[2] - res[1]^2
}

dph <- function(x, ph = ph(1), log = FALSE) {
  inv <- order(order(x))
  x <- c(0,sort(x))
  res <- mexp(t=x, transpose=TRUE, x=ph@alpha, A=ph@Q)
  y <- as.vector(ph@xi %*% res$result[,2:length(x)])[inv]
  if (log) {
    log(y)
  } else {
    y
  }
}

pph <- function(q, ph = ph(1), lower.tail = TRUE, log.p = FALSE) {
  inv <- order(order(q))
  x <- c(0,sort(q))
  if (lower.tail) {
    al <- c(ph@alpha,0)
    A <- diag.padding.zero(rbind2(cbind2(ph@Q,ph@xi),rep(0,ph@size+1)))
    res <- mexp(t=x, transpose=TRUE, x=al, A=A)
    y <- res$result[ph@size+1,2:length(x)][inv]
    if (log.p) {
      log(y)
    } else {
      y
    }
  } else {
    res <- mexp(t=x, transpose=TRUE, x=ph@alpha, A=ph@Q)
    y <- as.vector(rep(1,ph@size) %*% res$result[,2:length(x)])[inv]
    if (log.p) {
      log(y)
    } else {
      y
    }
  }
}

rph <- function(n, ph = ph(1)) {
  if (class(ph) == "cf1") {
    cf1.sample(n, ph)
  } else {
    sapply(1:n, function(a) ph.sample(ph))
  }
}

ph.sample <- function(ph) {
  s <- which(as.vector(rmultinom(n=1, size=1, prob=c(ph@alpha,0)))==1)
  t <- 0
  while (s != ph@size+1) {
    x <- c(ph@Q[s,], ph@xi[s])
    r <- -x[s]
    p <- x / r
    p[s] <- p[s] + 1
    t <- t + rexp(n=1, rate=r)
    s <- which(as.vector(rmultinom(n=1, size=1, prob=p))==1)
  }
  t
}

## S4 methods

setMethod("emfit.df", signature(model = "ph"),
  function(model, ...) {
    model@df
  }
)

setMethod("emfit.print", signature(model = "ph"),
  function(model, ...) {
    cat(gettextf("Size : %d\n", model@size))
    cat("Initial : ", model@alpha, "\n")
    cat("Exit    : ", model@xi, "\n")
    cat("Infinitesimal generator : \n")
    print(model@Q)
  }
)

setMethod("emfit.init", signature(model = "ph"),
  function(model, data, verbose = list(), ...) {
    ph.param.random(size=model@size, data=data,
      skelpi=model@alpha, skelQ=as.matrix(model@Q),
      skelxi=model@xi, verbose=verbose, class=class(model@Q))
  }
)

ph.param.random <- function(size, data, skelpi, skelQ, skelxi, verbose, class) {
  if (missing(size)) size <- length(skelpi)
  if (missing(skelpi)) skelpi <- rep(1, size)
  if (missing(skelQ)) skelQ <- matrix(1, size, size)
  if (missing(skelxi)) skelxi <- rep(1, size)

  mean <- mapfit.mean(data)

  alpha <- skelpi * runif(size)
  alpha <- alpha / sum(alpha)

  diag(skelQ) <- 0
  Q <- skelQ * matrix(runif(size*size), size, size)
  xi <- skelxi * runif(size)
  diag(Q) <- -(apply(Q, 1, sum) + xi)

  p <- ph(alpha=alpha, Q=Q, xi=xi, class=class)
  m <- ph.mean(p) / mean
  p@Q <- as(as.matrix(p@Q * m), class)
  p@xi <- p@xi * m
  p
}

#### estep

setMethod("emfit.estep", signature(model = "ph", data = "phdata.wtime"),
  function(model, data, ufact = 1.01, eps = sqrt(.Machine$double.eps), ...) {
    res <- .Call(phfit_estep_gen_wtime, model, data, eps, ufact)
    list(eres=list(etotal=res[[1]], eb=res[[2]], ey=res[[3]], ez=res[[4]], en=res[[5]]), llf=res[[6]])
  })

setMethod("emfit.estep", signature(model = "ph", data = "phdata.group"),
  function(model, data, ufact = 1.01, eps = sqrt(.Machine$double.eps), ...) {

  data@data$instant[is.na(data@data$counts)] <- 0
  data@data$counts[is.na(data@data$counts)] <- -1
  l <- data@size
  if (is.infinite(data@data$time[l])) {
    gdatlast <- data@data$counts[l]
    data@data <- data@data[-l,]
    data@size <- data@size - 1
  } else {
    gdatlast <- 0
  }

  ba <- msolve(alpha=1.0, A=-as.matrix(model@Q), x=model@alpha, transpose=TRUE)

  res <- .Call(phfit_estep_gen_group, model, ba, data, gdatlast, eps, ufact)
  list(eres=list(etotal=res[[1]], eb=res[[2]], ey=res[[3]], ez=res[[4]], en=res[[5]]), llf=res[[6]])
  })

#### mstep

setMethod("emfit.mstep", signature(model = "ph"),
  function(model, eres, data, ...) {
    res <- .Call(phfit_mstep_gen, model, eres, data)
    model@alpha <- res[[1]]
    model@xi <- res[[2]]
    model@Q@x <- res[[3]]
    model
  })

