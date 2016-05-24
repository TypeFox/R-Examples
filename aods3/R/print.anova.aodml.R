print.anova.aodml <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  mod <- x$models
  dfr <- x$anova.table
  nam <- rownames(dfr)
  cat("Likelihood ratio tests\n\n")
  sapply(mod, function(x) cat(x, "\n"))
  cat("\n")
## Modif RL 16/06/2013
##  List <- lapply(dfr, function(x) ifelse(is.na(x), "", format(x, digits = 4)))
##  dfr <- as.data.frame(t(do.call("rbind", List)))
  rownames(dfr) <- nam
## Ajout argument digits
  print(dfr, digits = digits)
  invisible(dfr)
  }
