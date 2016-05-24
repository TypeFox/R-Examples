### accumulate values from piecewise constant function of time
const2lin <- function(const, time = "time", stratum = NULL)
{
  if (is.character(stratum)) stratum <- match(stratum, colnames(const))
  if (is.character(time)) time <- match(time, colnames(const))
  m <- max(1, ncol(const[, -c(time, stratum)]))
  f <- function(x) c(0, -x[-length(x)]) + x
  g <- function(x) if (length(dim(x))) apply(x, 2, f) else f(x)
  h <- function(x) if (length(dim(x))) apply(x, 2, cumsum) else cumsum(x)
  if (!is.null(stratum)) {
    len <- by(const[, time], const[, stratum], g, simplify = FALSE)
    len <- do.call("rbind", lapply(len, matrix, ncol = m))
    lin <- by(const[, -c(time, stratum)] * len, const[, stratum], h,
              simplify = FALSE)
    lin <- do.call("rbind", lapply(lin, matrix, ncol = m))
  }
  else {
    len <- g(const[, time])
    lin <- h(const[, -time] * len)
  }
  lin <- cbind(lin, const[, c(time, stratum)])
  rownames(lin) <- rownames(const)
  colnames(lin) <- colnames(const)
  if (is.data.frame(const)) data.frame(lin)
  else lin
}
