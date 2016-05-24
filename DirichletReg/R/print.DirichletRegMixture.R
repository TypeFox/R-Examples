print.DirichletRegMixture <- function(x, ...){

  cat("\n\nDirichlet Mixture of",x$classes,"Latent Classes\n\n\nGroup Sizes and Probabilities:\n")

  sp.tab <- as.table(rbind(format(round(x$groups,0)),
                           format(round(x$groups/sum(x$groups),3))))
  dimnames(sp.tab) <- list(c("Size","Probability"), paste("Class", 1:x$classes))
  print(sp.tab, print.gap=2, justify="right")

  cat("\n\nLog-Likelihood:", round(x$logLik,3), "with", x$npar, "Parameters")
  cat("\nBIC:", round(BIC(x),1))
  cat("\nEntropy:", round(x$entropy,3),"\n\n")
}
