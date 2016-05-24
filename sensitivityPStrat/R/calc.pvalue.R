calc.pvalue <- function(x, var, test=c(upper=FALSE,lower=FALSE,twoSided=TRUE)) {
  a <- pnorm(x/sqrt(var))

  cbind(if(test['upper']) 1 - a,
        if(test['lower']) a,
        if(test['twoSided']) 2 * ifelse(a > 0.5, 1 - a, a))
}
