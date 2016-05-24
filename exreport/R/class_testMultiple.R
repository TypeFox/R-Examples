
is.testMultiple <- function(x) {
  is(x, "testMultiple")
}

print.testMultiple <- function(x, ...) {
  print(x$f)
  cat(sprintf("%s post-hoc test with %s p-value adjustment for output %s", x$tags$scope, x$tags$method, x$tags$target))
}

.testMultiple <- function(names, pvalues, friedman, experiment, alpha, tags){
  
  ph <- list(
    "names"       = names,
    "pvalues"     = pvalues,
    "friedman"    = friedman,
    "experiment"  = experiment,
    "tags"        = tags
  )
  
  class(ph) <- c("testMultiple", "reportable")
  ph
}



