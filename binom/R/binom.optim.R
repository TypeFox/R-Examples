binom.optim <- function(n, conf.level = 0.95, method = binom.lrt, k = n%/%2 + 1, p0 = 0, transform = TRUE,
                        plot = FALSE, tol = .Machine$double.eps^0.5, start = NULL, ...) {
  coverage <- conf.level
  if(is.character(method) && exists(method)) method <- get(method)
  if((k < 1 || k > n%/%2 + 1))
    stop(sprintf("k must be between 1 and %d", n%/%2 + 1))
  link <- if(transform) binomial()$linkfun else I
  linkinv <- if(transform) binomial()$linkinv else I
  lower <- link(rep(tol, k))
  upper <- link(rep(0.10, k))
  alpha <- if(is.null(start)) {
    link(rep(1 - conf.level, length.out = k))
  } else {
    link(rep(1 - start, length.out = k))
  }
  fit <- optim(alpha, binom.optim.objective, method = "L-BFGS-B",
               n = n, k = k, p0 = p0, ci.method = method,
               coverage = conf.level, transform = linkinv, ...)
  fit$par <- 1 - linkinv(fit$par)
  conf.level <- c(fit$par, rep(conf.level, n%/%2 + 1 - k))
  conf.level <- c(conf.level, rev(conf.level))
  if(n%%2 == 0) conf.level <- conf.level[-(n%/%2 + 1)]
  if(plot) {
    print(binom.plot(n, method, conf.level = conf.level, actual = coverage, type = "xy", ...))
  }
  fit$confint <- cbind(method(0:n, n, conf.level, ...), conf.level = conf.level)
  fit
}

binom.optim.objective <- function(alpha, n, k, p0, ci.method, coverage, transform, half = TRUE, ...) {
  alpha <- transform(alpha)
  conf.level <- 1 - c(alpha, rep(1 - coverage, n%/%2 + 1 - k))
  k <- length(conf.level)
  if(k != 1 && k < n + 1) {
    conf.level <- c(conf.level, rev(conf.level))
    if(n%%2 == 0) conf.level <- conf.level[-k]
  }
  ci <- ci.method(0:n, n, conf.level = conf.level, ...)[c("x", "lower", "upper")]
  b <- sort(unique(unlist(ci[-1])))
  if(half) b <- c(b[b < 0.5], 0.5)
  mid.b <- b[-length(b)] + diff(b) * 0.5
  v <- 0
  for(i in seq(along = mid.b)) {
    x <- ci$x[mid.b[i] > ci$lower & mid.b[i] < ci$upper]
    if(!length(x))
      stop("optimizer is divergent. Use a different k or a different method.")
    b0 <- max(p0, b[i])
    b1 <- min(1 - p0, b[i + 1])
    if(b1 < b0) next
#    print(c(b0 = b0, b1 = b1))
    v <- v + integrate.poly(x, n, b0, b1, coverage, "pbeta")
  }
  v * if(half) 2 else 1
}

integrate.poly <- function(x, n, lower = 0, upper = 1, conf.level = 0.95, method = c("pbeta", "polynom")) {
  method <- match.arg(method)
  if(length(n) > 1) stop("n must be at most length 1")
  if(length(x) < 1) stop("x must be at least length 1")
  bounds <- sort(c(lower[1], upper[1]))
  x <- sort(x)
  if(method == "polynom") {
    require(polynom)
    integrand <- polynomial(0)
    p1 <- polynomial(c(0, 1))
    p2 <- polynomial(c(1, -1))
    for(k in x) integrand <- integrand + p1^k * p2^(n - k) * choose(n, k)
    diff(predict(integral((integrand - conf.level)^2), bounds))
  } else {
    pbeta2 <- function(p, a, b, r) {
      pr <- exp(r + sapply(p, pbeta, shape1 = a, shape2 = b, log.p = TRUE))
      if(is.matrix(pr)) pr[, 2] - pr[, 1] else pr[2] - pr[1]
    }
    m <- length(x)
    log.np1 <- log(n + 1)
    ## squared terms: (choose(n, x) * p^x * (1 - p)^(n - x))^2
    a <- 2 * x + 1
    b <- 2 * (n - x) + 1
    rx <- lbeta(x + 1, n - x + 1)
    r <- lbeta(a, b) - 2 * (rx + log.np1)
    pr.1 <- sum(pbeta2(bounds, a, b, r))
    ## cross terms: -2 * (1 - alpha) * choose(n, x) * p^x * (1 - p)^(n - x)
    a <- x + 1
    b <- n - x + 1
    r <- -log.np1
    pr.2 <- sum(-2 * conf.level * pbeta2(bounds, a, b, r))
    ## cross terms: 2 * choose(n, x) * p^x * (1 - p)^(n - x) * choose(n, y) * p^y * (1 - p)^(n - y)
    pr.3 <- if(m > 1) {
      y <- unlist(lapply(seq(m - 1), function(i) x[-seq(i)]))
      x <- rep(x[-m], rev(seq(m - 1)))
      a <- x + y + 1
      b <- 2 * n - x - y + 1
      rx <- rep(rx[-m], rev(seq(length(rx) - 1)))
      r <- lbeta(a, b) - (rx + lbeta(y + 1, n - y + 1) + 2 * log.np1)
      sum(2 * pbeta2(bounds, a, b, r))
    } else {
      0
    }
    ## Un-vectorized version (slow for large length(x))
    ##    pr.1 <- pr.2 <- pr.3 <- 0
    ##    log.np1 <- log(n + 1)
    ##    index.i <- seq(along = x)
    ##    for(i in index.i) {
    ##      x.i <- x[i]
    ##      a <- 2 * x.i + 1
    ##      b <- 2 * (n - x.i) + 1
    ##      r.i <- lbeta(x.i + 1, n - x.i + 1)
    ##      r <- lbeta(a, b) - 2 * (r.i + log.np1)
    ##      pr.1 <- pr.1 + diff(exp(r + pbeta(bounds, a, b, log.p = TRUE)))
    ##      a <- x.i + 1
    ##      b <- n - x.i + 1
    ##      pr.2 <- pr.2 - 2 * conf.level * diff(exp(pbeta(bounds, a, b, log.p = TRUE) - log.np1))
    ##      index.j <- seq(along = x)[-1:-i]
    ##      for(j in index.j) {
    ##        x.j <- x[j]
    ##        a <- x.i + x.j + 1
    ##        b <- 2 * n - x.i - x.j + 1
    ##        r <- lbeta(a, b) - (r.i + lbeta(x.j + 1, n - x.j + 1) + 2 * log.np1)
    ##        pr.3 <- pr.3 + 2 * diff(exp(r + pbeta(bounds, a, b, log.p = TRUE)))
    ##      }
    ##    }
    pr.1 + pr.2 + pr.3 + diff(bounds) * conf.level^2
  }
}
