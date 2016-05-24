# 03/Mar/2011 Bug fix: pcX(nrow(M)) instead of pcX(nrow(M))[,-1] in thurstone

thurstone <- function(M){
  # Fit Thurstone-Mosteller model (Case V) using glm
  # M: paired-comparison matrix (absolute frequencies)
  # author: Florian Wickelmaier (wickelmaier@web.de)
  # last mod: 03/Mar/2011

  y1 <- t(M)[lower.tri(M)]
  y0 <- M[lower.tri(M)]

  tm.glm <- glm(cbind(y1, y0) ~ pcX(nrow(M)) - 1, binomial(probit))
  estimate <- c(0, coef(tm.glm))
  names(estimate) <- colnames(M)

  gof <- c("-2logL" = deviance(tm.glm), df = tm.glm$df.residual,
    pval = 1 - pchisq(deviance(tm.glm), tm.glm$df.residual))

  z <- list(estimate=estimate, goodness.of.fit=gof, tm.glm=tm.glm)
  class(z) <- "thurstone"
  z
}


print.thurstone <- function(x, digits=max(3, getOption("digits") - 3),
  na.print="", ...){
  cat("\nThurstone-Mosteller model (Case V)\n\n")
  cat("Parameter estimates:\n")
  print.default(format(x$estimate, digits = digits), print.gap = 2,
      quote = FALSE)
  G2   <- x$goodness.of.fit[1]
  df   <- x$goodness.of.fit[2]
  pval <- x$goodness.of.fit[3]
  cat("\nGoodness of fit (-2 log likelihood ratio):\n")
  cat("\tG2(", df, ") = ", format(G2, digits=digits), ", p = ",
      format(pval,digits=digits), "\n", sep="")
  cat("\n")
  invisible(x)
}
