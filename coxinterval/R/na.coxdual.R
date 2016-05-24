na.coxdual <- function(object, ...) UseMethod("na.coxdual")
na.coxdual.default <- function(object, ...)
{
  n <- length(object)
  omit <- FALSE
  vars <- seq_len(n)
  for(j in vars) {
    x <- object[[j]]
    if (!is.atomic(x))
      next
    if (inherits(x, "Surv")) {
      ## allow missing values for 'time1' or 'time2', but not both
      x[is.na(x[, 1]) & !is.na(x[, 2]), 1] <- -1
      x[!is.na(x[, 1]) & is.na(x[, 2]), 2] <- -1
    }
    ## save states attribute for 'trans' term
    ## allow missing values for ? -> 2 transitions
    if (!is.null(states <- attr(x, "states"))) {
      types <- attr(x, "types")
      x[is.na(x[, 1]) & x[, 2] %in% states[3], 1] <- -Inf
    }
    x <- is.na(x)
    d <- dim(x)
    if (is.null(d) || length(d) != 2L)
      omit <- omit | x
    else
      for(ii in 1L:d[2L])
        omit <- omit | x[, ii]
  }
  xx <- object[!omit, , drop = FALSE]
  if (any(omit > 0L)) {
    temp <- seq(omit)[omit]
    names(temp) <- attr(object, "row.names")[omit]
    attr(temp, "class") <- "exclude"
    attr(xx, "na.action") <- temp
  }
  # save states attribute for 'trans' term
  if (!is.null(states)) {
    attr(xx, "states") <- states
    attr(xx, "types") <- types
  }
  xx
}
