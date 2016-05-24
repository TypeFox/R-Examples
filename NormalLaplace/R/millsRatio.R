millsR <- function(y, log = FALSE) {

  if (log) {
    millsR <- pnorm(y, lower.tail = FALSE, log.p = TRUE) -
              dnorm(y, log = TRUE)
  } else {
    millsR <- ifelse(y < 30, pnorm(y, lower.tail = FALSE)/dnorm(y),
                     exp(pnorm(y, lower.tail = FALSE, log.p = TRUE) -
                         dnorm(y, log = TRUE)))
  }

  return(millsR)
}
