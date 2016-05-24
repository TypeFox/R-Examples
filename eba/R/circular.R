# Circular triads (transitivity violations)
# May/10/2010: Bug fix, p-value for zeta is 1 - p-value of chi-square test
#              (reported by Wolfgang Ellermeier)

circular <- function(mat){
  mat <- as.matrix(mat)
  diag(mat) <- 0         # remove diagonal
  n <- ncol(mat)
  if(n < 8) warning("Chi-square approximation might be incorrect.",
    call. = FALSE)
  T <- n * (n - 1) * (2*n - 1)/12 - .5 * sum(colSums(mat)^2)
  if(n%%2) T.max <- n*(n^2 - 1)/24
  else     T.max <- n*(n^2 - 4)/24
  zeta <- 1 - T/T.max
  df <- n * (n - 1) * (n - 2) / (n - 4)^2
  chi2 <- 8/(n - 4) * (.25*choose(n, 3) - T + .5) + df
  z <- list(T=T, T.max=T.max, zeta=zeta, chi2=chi2, df=df, p=pchisq(chi2,df))
  class(z) <- "circular"
  z
}


print.circular <- function(x, digits = max(3,getOption("digits")-4), ...){
  cat("\nCircular triads (intransitive cycles)\n\n")
  # cat("T = ",x$T, ", T.max = ",x$T.max, ", zeta = ",x$zeta,
  #     "\nchi2 = ",x$chi2, ", df = ",x$df, ", p-value = ",x$p, "\n",sep="")
  cat("T = ",x$T, ", T.max = ",x$T.max, ", zeta = ",x$zeta,
      ", p-value = ",x$p, "\n",sep="")
  cat("alternative hypothesis: intransitivity is not by chance\n")
  cat("\n")
  invisible(x)
}
