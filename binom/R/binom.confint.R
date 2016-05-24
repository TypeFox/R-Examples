binom.methods <- c("wilson", "agresti-coull", "ac", "exact", "prop.test", "profile", "lrt",
                   "asymptotic", "bayes", "cloglog", "logit", "probit", "all")

binom.confint <- function(x, n, conf.level = 0.95, methods = "all", ...) {
  if((length(x) != length(n))) {
    m <- cbind(x = x, n = n)
    x <- m[, "x"]
    n <- m[, "n"]
  }
  res <- NULL
  method <- pmatch(methods, binom.methods)
  if(all(is.na(method))) {
    methods <- paste(paste("\"", methods, "\"", sep = ""), collapse = ", ")
    stop("method(s) ", methods, "\" not matched")
  }
  if(any(is.na(method))) {
    warning("method(s) ", methods[is.na(method)], " not matched")
    method <- methods[!is.na(method)]
  }
  method <- binom.methods[method]
  out <- x > n | x < 0
  if(any(out)) {
    warning("observations with more successes than trials detected and removed")
    x <- x[!out]
    n <- n[!out]
  }
  xn <- data.frame(x = x, n = n)
  all.methods <- any(method == "all")
  p <- x/n
  alpha <- 1 - conf.level
  alpha <- rep(alpha, length = length(p))
  alpha2 <- 0.5 * alpha
  z <- qnorm(1 - alpha2)
  z2 <- z * z
  res <- NULL
  if(any(method %in% c("agresti-coull", "ac")) || all.methods) {
    .x <- x + 0.5 * z2
    .n <- n + z2
    .p <- .x/.n
    lcl <- .p - z * sqrt(.p * (1 - .p)/.n)
    ucl <- .p + z * sqrt(.p * (1 - .p)/.n)
    res.ac <- data.frame(method = rep("agresti-coull", NROW(x)),
                         xn, mean = p, lower = lcl, upper = ucl)
    res <- res.ac
  }
  if(any(method == "asymptotic") || all.methods) {
    se <- sqrt(p * (1 - p)/n)
    lcl <- p - z * se
    ucl <- p + z * se
    res.asymp <- data.frame(method = rep("asymptotic", NROW(x)),
                            xn, mean = p, lower = lcl, upper = ucl)
    res <- if(is.null(res)) res.asymp else rbind(res, res.asymp)
  }
  if(any(method == "bayes") || all.methods) {
    res.bayes <- binom.bayes(x, n, conf.level = conf.level, ...)
    res.bayes <- res.bayes[c("method", "x", "n", "mean", "lower", "upper")]
    res <- if(is.null(res.bayes)) res.bayes else rbind(res, res.bayes)
  }
  if(any(method == "cloglog") || all.methods) {
    if(any(method != "exact")) {
      x1 <- x == 0
      x2 <- x == n
    }
    inner <- !x1 & !x2
    log.mu <- sd <- lcl <- ucl <- rep(NA, length(x))
    log.mu[inner] <- log(-log(p[inner]))
    sd[inner] <- sqrt(var.cloglog(p[inner], n[inner]))
    lcl[inner] <- exp(-exp(log.mu[inner] + z[inner] * sd[inner]))
    ucl[inner] <- exp(-exp(log.mu[inner] - z[inner] * sd[inner]))
    lcl[x1] <- rep(0, sum(x1))
    lcl[x2] <- alpha2[x2]^(1/n[x2])
    ucl[x1] <- 1 - alpha2[x1]^(1/n[x1])
    ucl[x2] <- rep(1, sum(x2))
    res.cloglog <- data.frame(method = rep("cloglog", NROW(x)),
                              xn, mean = p, lower = lcl, upper = ucl)
    res <- if(is.null(res)) res.cloglog else rbind(res, res.cloglog)
  }
  if(any(method == "exact") || all.methods) {
    x1 <- x == 0
    x2 <- x == n
    lb <- ub <- x
    lb[x1] <- 1
    ub[x2] <- n[x2] - 1
    lcl <- 1 - qbeta(1 - alpha2, n + 1 - x, lb)
    ucl <- 1 - qbeta(alpha2, n - ub, x + 1)
    if(any(x1)) lcl[x1] <- rep(0, sum(x1))
    if(any(x2)) ucl[x2] <- rep(1, sum(x2))
    res.exact <- data.frame(method = rep("exact", NROW(x)),
                            xn, mean = p, lower = lcl, upper = ucl)
    res <- if(is.null(res)) res.exact else rbind(res, res.exact)
  }
  if(any(method == "logit") || all.methods) {
    if(any(method != "exact")) {
      x1 <- x == 0
      x2 <- x == n
    }
    inner <- !x1 & !x2
    logit.p <- sd <- lcl <- ucl <- rep(NA, length(x))
    logit.p[inner] <- log(p[inner]) - log1p(-p[inner])
    sd[inner] <- sqrt(var.logit(p[inner], n[inner]))
    .lcl <- exp(logit.p[inner] - z[inner] * sd[inner])
    .ucl <- exp(logit.p[inner] + z[inner] * sd[inner])
    lcl[inner] <- .lcl/(1 + .lcl)
    ucl[inner] <- .ucl/(1 + .ucl)
    lcl[x1] <- rep(0, sum(x1))
    lcl[x2] <- alpha2[x2]^(1/n[x2])
    ucl[x1] <- 1 - alpha2[x1]^(1/n[x1])
    ucl[x2] <- rep(1, sum(x2))
    res.logit <- data.frame(method = rep("logit", NROW(x)),
                            xn, mean = p, lower = lcl, upper = ucl)
    res <- if(is.null(res)) res.logit else rbind(res, res.logit)
  }
  if(any(method == "probit") || all.methods) {
    if(any(method != "exact")) {
      x1 <- x == 0
      x2 <- x == n
    }
    inner <- !x1 & !x2
    probit.p <- sd <- lcl <- ucl <- rep(NA, length(x))
    probit.p[inner] <- qnorm(p[inner])
    sd[inner] <- sqrt(var.probit(p[inner], n[inner]))
    lcl[inner] <- pnorm(probit.p[inner] - z[inner] * sd[inner])
    ucl[inner] <- pnorm(probit.p[inner] + z[inner] * sd[inner])
    lcl[x1] <- rep(0, sum(x1))
    lcl[x2] <- alpha2[x2]^(1/n[x2])
    ucl[x1] <- 1 - alpha2[x1]^(1/n[x1])
    ucl[x2] <- rep(1, sum(x2))
    res.probit <- data.frame(method = rep("probit", NROW(x)),
                             xn, mean = p, lower = lcl, upper = ucl)
    res <- if(is.null(res)) res.probit else rbind(res, res.probit)
  }
  if(any(method == "profile") || all.methods) {
    res.prof <- binom.profile(x, n, conf.level, ...)
    res <- if(is.null(res)) res.prof else rbind(res, res.prof)
  }
  if(any(method == "lrt") || all.methods) {
    res.lrt <- binom.lrt(x, n, conf.level = conf.level, ...)
    res <- if(is.null(res)) res.lrt else rbind(res, res.lrt)
  }
  if(any(method == "prop.test") || all.methods) {
    ci <- lapply(seq_along(x), function(i) stats::prop.test(x[i], n[i])$conf.int)
    lcl <- sapply(ci, "[", 1)
    ucl <- sapply(ci, "[", 2)
    res.prop.test <- data.frame(method = rep("prop.test", NROW(x)),
                                xn, mean = p, lower = lcl, upper = ucl)
    res <- if(is.null(res)) res.prop.test else rbind(res, res.prop.test)
  }
  if(any(method == "wilson") || all.methods) {
    p1 <- p + 0.5 * z2/n
    p2 <- z * sqrt((p * (1 - p) + 0.25 * z2/n)/n)
    p3 <- 1 + z2/n
    lcl <- (p1 - p2)/p3
    ucl <- (p1 + p2)/p3
    # x1 <- x == 1
    # x2 <- x == n - 1
    # if(any(x1)) lcl[x1] <- -log(1 - alpha[x1])/n[x1]
    # if(any(x2)) ucl[x2] <- 1 + log(1 - alpha[x2])/n[x2]
    res.wilson <- cbind(method = rep("wilson", NROW(x)),
                        xn, mean = p, lower = lcl, upper = ucl)
    res <- if(is.null(res)) res.wilson else rbind(res, res.wilson)
  }
  attr(res, "conf.level") <- conf.level
  row.names(res) <- seq(nrow(res))
  res
}

binom.wilson <- function(x, n, conf.level = 0.95, ...)
  binom.confint(x, n, conf.level, "wilson")

binom.agresti.coull <- function(x, n, conf.level = 0.95, ...)
  binom.confint(x, n, conf.level, "ac")

binom.exact <- function(x, n, conf.level = 0.95, ...)
  binom.confint(x, n, conf.level, "exact")

binom.asymp <- function(x, n, conf.level = 0.95, ...)
  binom.confint(x, n, conf.level, "asymp")

binom.prop.test <- function(x, n, conf.level = 0.95, ...)
  binom.confint(x, n, conf.level, "prop.test")

binom.cloglog <- function(x, n, conf.level = 0.95, ...)
  binom.confint(x, n, conf.level, "cloglog")

binom.logit <- function(x, n, conf.level = 0.95, ...)
  binom.confint(x, n, conf.level, "logit")

binom.probit <- function(x, n, conf.level = 0.95, ...)
  binom.confint(x, n, conf.level, "probit")
