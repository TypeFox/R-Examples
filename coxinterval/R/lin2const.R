### convert piecewise linear values to piecewise constant values
lin2const <- function(lin, time = "time", stratum = NULL)
{
  if (is.character(stratum)) stratum <- match(stratum, colnames(lin))
  if (is.character(time)) time <- match(time, colnames(lin))
  m <- max(1, ncol(lin[, -c(time, stratum)]))
  f <- function(x) c(0, -x[-length(x)]) + x
  g <- function(x) if (length(dim(x))) apply(x, 2, f) else f(x)
  if (!is.null(stratum)) {
    jump <- by(lin[, -c(time, stratum)], lin[, stratum], g, simplify = FALSE)
    jump <- do.call(rbind, lapply(jump, matrix, ncol = m))
    len <- by(lin[, time], lin[, stratum], g, simplify = FALSE)
    len <- do.call("rbind", lapply(len, matrix, ncol = m))
  }
  else {
    jump <- g(lin[, -time])
    len <- g(lin[, time])
  }
  const <- jump
  const[len > 0] <- const[len > 0] / len[len > 0]
  const <- cbind(const, lin[, c(time, stratum)])
  rownames(const) <- rownames(lin)
  colnames(const) <- colnames(lin)
  if (is.data.frame(lin)) data.frame(const)
  else const
}
