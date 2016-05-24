### convert jumps to cumulative values
jump2step <- function(jump, time = "time", stratum = NULL)
{
  if (is.character(stratum)) stratum <- match(stratum, colnames(jump))
  if (is.character(time)) time <- match(time, colnames(jump))
  m <- max(1, ncol(jump[, -c(time, stratum)]))
  f <- function(x) if (length(dim(x))) apply(x, 2, cumsum) else cumsum(x)
  if (!is.null(stratum)) {
    step <- by(jump[, -c(time, stratum)], jump[, stratum], f, simplify = FALSE)
    step <- do.call("rbind", lapply(step, matrix, ncol = m))
  }
  else step <- f(jump[, -time])
  step <- cbind(step, jump[, c(time, stratum)])
  rownames(step) <- rownames(jump)
  colnames(step) <- colnames(jump)
  if (is.data.frame(jump)) data.frame(step)
  else step
}
