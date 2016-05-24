print.graf <-
function(x, ...){
  # print function for graf objects
  cat("GRaF model fitted with\n")
  cat(paste("  ", sum(x$obsy), "presence records\n"))
  cat(paste("  ", sum(1 - x$obsy), "absence records"))
  cat(paste("\n   and", ncol(x$obsx) , "covariates\n"))
}
