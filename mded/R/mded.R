mded <-
function(distr1, distr2, detail = FALSE, independent = TRUE)
{
  distr.names <- c(deparse(substitute(distr1)),
                  deparse(substitute(distr2)))

  if(independent == TRUE) {
    n1 <- length(distr1)
    n2 <- length(distr2)
    tn <- n1 * n2
    if(detail == FALSE) {
      cond <- 0
      for (i in 1:n2) {
        cond <- cond + sum((distr1 - distr2[i]) <= 0)
      }
    } else {
      dif <- rep(distr1, times = n2) - rep(distr2, each = n1)
      cond <- sum(dif <= 0)
    }
  } else {
    stopifnot(length(distr1) == length(distr2))
    tn <- length(distr1)
    dif <- distr1 - distr2
    cond <- sum(dif <= 0)
  }
 
  out <- list(stat   = c(cond / tn),
              means  = c(mean(distr1), mean(distr2)),
              cases  = c(cond, tn - cond),
              distr1 = distr1,
              distr2 = distr2,
              distr.names = distr.names
             )

  if(detail == TRUE) {
    out <- c(out, diff = list(dif))
  }

  class(out) <- "mded"

  return(out)
}
