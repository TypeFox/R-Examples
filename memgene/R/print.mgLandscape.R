
print.mgLandscape <- function(x, ...) {

  interpretation <- data.frame(model=names(x$summary)[c(2,4,6,8,9)], description=c(
                                  "explained by spatial predictors",
                                  "spatial and explained by selected patterns in the model",
                                  "spatial and explained by coordinates not patterns in the model",
                                  "spatial and confounded between the model and coordinates",
                                  "residual and not explained by spatial predictors"), stringsAsFactors=FALSE)
  interpretation <- apply(format(interpretation, justify="left"), 1, function(y) paste0(y[1], "\t", y[2], "\n"))
  cat("mgLandscape Analysis\n")
  print(x$summary[, 1:9], digits=3, row.names=FALSE, ...)
  cat("\n")
  cat("Interpretation:\n")
  cat("Proportion of variation in genetic distance that is... (RsqAdj)\n")
  sapply(interpretation, function(x) cat(x))
  cat("\n")

}