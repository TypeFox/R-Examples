
seasonal.dummies <- function(x)
{
  S <- frequency(x)
  SD <- do.call("rbind", replicate(ceiling(length(x)/S), diag(S), simplify = FALSE))
  SD <- ts(SD, frequency = S, start = c(start(x)[1], 1))
  # ignore warning about 'start' or 'end' "value not changed"
  SD <- suppressWarnings(window(SD, start = start(x), end = end(x)))
  colnames(SD) <- paste("SD", seq_len(S), sep = "")
  SD
}

seasonal.cycles <- function(x)
{
  n <- length(x)
  S <- frequency(x)
  tmp <- matrix(seq_len(n), nrow = floor(0.5 * S - 1), ncol = n, byrow = TRUE)
  seqsm1 <- seq_len(nrow(tmp))
  tmp <- (2 * seqsm1 * pi / S) * tmp
  SD <- rbind(cos(tmp), sin(tmp))
  SD <- t(SD[c(rbind(seqsm1, seqsm1 + nrow(tmp))),])
  if ((S %% 2) == 0)
    SD <- cbind(SD, rep(c(-1, 1), len = n))
  SD <- ts(SD, frequency = S, start = c(start(x)[1], 1))
  colnames(SD) <- paste("SD", seq_len(ncol(SD)), sep = "")
  SD
}
