propCI <-
function(x, n, method = "all", level = 0.95, sortby = "level"){
  ## Check 'x' and 'n'
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (length(x) != 1 && length(n) != 1 && length(x) != length(n))
    stop("'x' and 'n' cannot have different lengths")

  if (length(x) > length(n)) n <- rep(n, length(x))
  if (length(n) > length(x)) x <- rep(x, length(n))

  ## Define list of methods
  m <- character()
  if (any(method == "all")){
    method <- c("agresti.coull", "exact", "jeffreys", "wald", "wilson")
    m <- method
  } else {
    for (i in seq_along(method)){
      if (any(method[i] == c("agresti.coull", "agresti-coull", "ac")))
        m[i] <- "agresti.coull"
      if (any(method[i] == c("asymptotic", "normal", "wald")))
        m[i] <- "wald"
      if (any(method[i] == c("clopper-pearson", "cp", "exact")))
        m[i] <- "exact"
      if (any(method[i] == c("jeffreys", "bayes")))
        m[i] <- "jeffreys"
      if (any(method[i] == c("wilson")))
        m[i] <- "wilson"
    }
  }

  if (length(m) == 0 | any(m == ""))
    stop(paste("'method' must be",
               "agresti.coull, exact, jeffreys, wald, or wilson"))

  ## Check 'level'
  checkInput(level, "level", range = c(0, 1))

  ## Check 'sortby'
  checkInput(sortby, "sortby", value = c("level", "method"))

  ## Create data.frame
  if (sortby == "level"){
    prm <- level
    sec <- method
    m <- rep(m, times = length(x) * length(level))
    mm <- rep(method, times = length(x) * length(level))
    p <- rep(level, each = length(x) * length(method))
    x <- rep(x, each = length(method), times = length(level))
    n <- rep(n, each = length(method), times = length(level))
  } else {
    prm <- method
    sec <- level
    m <- rep(m, each = length(x) * length(level))
    mm <- rep(method, each = length(x) * length(level))
    p <- rep(level, times = length(x) * length(method))
    x <- rep(x, each = length(level), times = length(method))
    n <- rep(n, each = length(level), times = length(method))
  }

  data <- data.frame(x = x, n = n, q = x/n, m = mm, p = p, l = NA, u = NA)

  ## Estimate intervals
  for (i in seq(nrow(data))){
    l <- ciLevel(data$x[i], data$n[i], data$p[i])
    cl <- switch(EXPR = as.character(m[i]),
      agresti.coull  = propCI_agresticoull(data$x[i], data$n[i], l),
      exact          = propCI_exact(data$x[i], data$n[i], l),
      jeffreys       = propCI_jeffreys(data$x[i], data$n[i], l),
      wald           = propCI_wald(data$x[i], data$n[i], l),
      wilson         = propCI_wilson(data$x[i], data$n[i], l))
    data$l[i] <- cl[1]
    data$u[i] <- cl[2]
  }

  names(data) <- c("x", "n", "p", "method", "level", "lower", "upper")
  return(data)
}
