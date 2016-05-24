
is.testMultiplePairwise <- function(x) {
  is(x, "testMultiplePairwise")
}

#' @export
print.testMultiplePairwise <- function (x, ...) {
  print.testMultiple(x)  
}

#' @export
summary.testMultiplePairwise <- function (object, ...) {
  cat("---------------------------------------------------------------------\n")
  summary(object$friedman)
  cat("---------------------------------------------------------------------\n")
  cat(sprintf("Pairwise post hoc test for output %s\n", object$tags$target))
  cat(sprintf("Adjust method: %s\n", object$tags$method))
  cat(sprintf("alpha = %.4f\n", object$alpha))
  cat("\n")
  cat("p-values:\n")
  d <- cbind(object$names, p=object$pvalue)
  means <- rowMeans(object$friedman$ranks)
  means <- means[order(means)]
  t <- reshape2::dcast(d, method1~method2, value.var="p")
  me <- t[,1]
  t <- t[,-1]
  rownames(t) <- me
  print(t[names(means)[-length(means)],names(means)[-1]])  
  cat("---------------------------------------------------------------------\n")
}

#Anonymous constructor
.testMultiplePairwise <- function(names, pvalues, friedman, experiment, alpha, target, scope, method, tags) {

  newTags <- .metaTags(alpha   = alpha,
                       target  = target,
                       scope   = scope,
                       method  = method,
                       alias   = "TestMultiplePairwise")
  tags <- .updateTags(tags, newTags)
  
  #Build the superclass object
  ph <- .testMultiple(names      = names, 
                      pvalues    = pvalues, 
                      friedman   = friedman,
                      experiment = experiment, 
                      tags       = tags
                      )
  
  #Specifical field for this class
  class(ph) <- append("testMultiplePairwise",class(ph))
  
  ph
}
