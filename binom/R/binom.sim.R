binom.sim <- function(M = 200, n = 100, p = 0.5,
                      conf.level = 0.95, methods = "all", ...) {
  .captured <- function(x, p)
    c(x = sum(p >= x$lower & p <= x$upper), width = mean(x$width))
  x <- rbinom(M, n, p)
  ci <- binom.confint(x, n, conf.level, methods, ...)
  ci$width <- with(ci, upper - lower)
  cap <- lapply(split(ci, ci$method), .captured, p = p)
  res <- list()
  for(method in names(cap)) {
    res[[method]] <- binom.confint(cap[[method]]["x"], M, conf.level, method)
    res[[method]] <- cbind(res[[method]], width = cap[[method]]["width"])
  }
  res <- do.call("rbind", res)[c("mean", "lower", "upper", "width")]
  attr(res, "call") <- match.call()
  res
}
