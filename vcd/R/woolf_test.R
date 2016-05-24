woolf_test <- function(x) {
  DNAME <- deparse(substitute(x))
  if (any(x == 0))
    x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1] * x[2,2]) / (x[1,2] * x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  o <- log(or)
  e <- weighted.mean(log(or), w)
  STATISTIC <- sum(w * (o - e)^2)
  PARAMETER <- k - 1
  PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
  METHOD <- "Woolf-test on Homogeneity of Odds Ratios (no 3-Way assoc.)"
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  structure(list(statistic = STATISTIC, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name = DNAME, observed = o, 
                 expected = e), class = "htest")
}

