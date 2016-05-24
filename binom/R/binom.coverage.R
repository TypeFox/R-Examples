binom.coverage <- function(p, n, conf.level = 0.95, method = "all", ...) {
  if(missing(p)) p <- seq(0, 1, length = 200)
  x <- unlist(lapply(lapply(n, ":", 0), rev))
  n <- rep(n, n + 1)
  ci <- if(is.function(method)) {
    method(x, n, conf.level = conf.level, ...)
  } else if (method %in% binom.methods) {
    binom.confint(x, n, conf.level, method, ...)
  } else if(is.character(method) && exists(method) && method != "all") {
    get(method)(x, n, conf.level = conf.level, ...)
  } else {
    stop("Don't know what to do with method '", method, "'.")
  }
  required.names <- c("method", "x", "n", "lower", "upper")
  if (!all(required.names %in% names(ci))) {
    msg <- "The following names are required in the return data for method '%s::%s':"
    msg <- paste0(msg, "\n", paste(required.names, collapse = ", "))
    method.name <- if (is.function(method)) deparse(substitute(method)) else method
    package.name <- sub("package:", "", find(method.name)[1L], fixed = TRUE)
    msg <- sprintf(msg, package.name, method.name)
    stop(msg)
  }
  ci <- ci[required.names]
  z <- merge(ci, data.frame(p = p))
  z$coverage <- with(z, (p >= lower & p <= upper) * dbinom(x, n, p))
  z <- aggregate(z["coverage"], z[c("method", "p", "n")], sum)
  z <- z[order(z$method, z$p), ]
  row.names(z) <- seq(NROW(z))
  z
}
