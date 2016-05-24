### convert cumulative values to a step function
step2stepfun <- function(step, time = "time", stratum = NULL) {
  if (is.character(stratum)) stratum <- match(stratum, colnames(step))
  if (is.character(time)) time <- match(time, colnames(step))
  f <- function(x, y) stepfun(y[-1], x)
  g <- function(x) {
    if (length(dim(x[, -time]))) apply(x[, -time], 2, f, y = x[, time])
    else f(x[, -time], x[, time])
  }
  if (!is.null(stratum))
    fun <- unlist(by(step[, -stratum], step[, stratum], g, simplify = FALSE))
  else fun <- g(step)
  fun
}
