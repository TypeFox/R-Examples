
is.testMultipleControl <- function(x) {
  is(x, "testMultipleControl")
}

#' @export
print.testMultipleControl <- function (x, ...) {
  print.testMultiple(x)  
}

#' @export
summary.testMultipleControl <- function (object, ...) {
  cat("---------------------------------------------------------------------\n")
  summary(object$friedman)
  cat("---------------------------------------------------------------------\n")
  cat(sprintf("Control post hoc test for output %s\n", object$tags$target))
  cat(sprintf("Multiple %s tests\n", object$tags$method2))
  cat(sprintf("Adjust method: %s\n", object$tags$method))
  cat(sprintf("alpha = %.4f\n", object$alpha))
  cat("\n")
  cat(sprintf("Control method: %s\n", object$control))
  cat("p-values:\n")
  d <- cbind(m=object$names,d=object$pvalues)
  d <- d[with(d, order(-d)),]
  for (r in 1:nrow(d)-1)
    cat(sprintf("%15s\t%.4f\n", d[r,1], d[r,2]))
  cat("---------------------------------------------------------------------\n")
}

#Anonymous constructor
.testMultipleControl <- function(names, control, pvalues, friedman, experiment,  alpha, target, scope, method, tags) {
  
  newTags <- .metaTags(alpha   = alpha,
                       target  = target,
                       scope   = scope,
                       method  = method,
                       alias   = "TestMultipleControl")
  tags <- .updateTags(tags, newTags)
  
  
  #Build the superclass object
  ph <- .testMultiple(names      = names, 
                      pvalues    = pvalues, 
                      friedman   = friedman,
                      experiment = experiment, 
                      tags       = tags
                      )
  
  # Specifical field for this class
  ph$control <- control
  
  class(ph)  <- append("testMultipleControl", class(ph))
  
  ph
}