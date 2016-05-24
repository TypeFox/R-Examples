vec2tfl <- function (x)
{
  freqs <- table(as.factor(x))
  idx <- order(freqs, decreasing=TRUE)
  result <- tfl(f=as.vector(freqs)[idx], type=names(freqs)[idx])
  if (N(result) != length(x))
    stop("internal error (inconsistency between sample size and type frequency list")
  result
}
