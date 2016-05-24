### convert cumulative values to jumps
step2jump <- function(step, time = "time", stratum = NULL)
{
  if (is.character(stratum)) stratum <- match(stratum, colnames(step))
  if (is.character(time)) time <- match(time, colnames(step))
  m <- max(1, ncol(step[, -c(time, stratum)]))
  f <- function(x) c(0, -x[-length(x)]) + x
  g <- function(x) if (length(dim(x))) apply(x, 2, f) else f(x)
  if (!is.null(stratum)) {
    jump <- by(step[, -c(time, stratum)], step[, stratum], g, simplify = FALSE)
    jump <- do.call("rbind", lapply(jump, matrix, ncol = m))
  }
  else jump <- g(step[, -time])
  jump <- cbind(jump, step[, c(time, stratum)])
  rownames(jump) <- rownames(step)
  colnames(jump) <- colnames(step)
  if (is.data.frame(step)) data.frame(jump)
  else jump
}
