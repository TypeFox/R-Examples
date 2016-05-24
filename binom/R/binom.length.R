binom.length <- function(p, n, conf.level = 0.95, method = "all", ...) {
  expect.length <- function(ci, x, n, p)
    sum((ci$upper - ci$lower) * dbinom(x, n, p))
  params <- expand.grid(p = p, n = n)
  len <- vector("list", NROW(params))
  for(i in seq(NROW(params))) {
    p.i <- params$p[i]
    n.i <- params$n[i]
    x <- 0:n.i
    if(is.function(method)) {
      ci <- method(x, n.i, conf.level, ...)
    } else {
      ci <- binom.confint(x, n.i, conf.level, method, ...)
    }
    ci <- split(ci[c("lower", "upper")], ci$method)
    expect <- lapply(ci, expect.length, x, n.i, p.i)
    len[[i]] <- data.frame(method = names(expect))
    len[[i]]$n <- n.i
    len[[i]]$p <- p.i
    len[[i]]$length <- unlist(expect)
  }
  len <- do.call("rbind", len)
  row.names(len) <- seq(NROW(len))
  len
}
