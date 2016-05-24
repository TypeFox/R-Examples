is.testPaired <- function(x) {
  is(x, "testPaired")
}

print.testPaired <- function (x, ...) {
  cat( sprintf("Wilcoxon test p-value (%s the output named %s), comparing %s versus %s: %.4e", x$tags$objetive, x$tags$target, x$bestMethod, x$worstMethod, x$tags$pvalue) )
}

summary.testPaired <- function (object, ...) {
  cat( sprintf("Wilcoxon test (%s the output named %s)\n\n", object$tags$objetive,object$tags$target) )
  cat( sprintf("Statistic distributed according to Wilcoxon Signed Rank: %.4f\n", object$statistic) )
  if(object$tags$pvalue<object$tags$alpha)
    cat( sprintf("%s vs %s. Test rejected: p-value: %.4e < %.4f", object$bestMethod, object$worstMethod, object$tags$pvalue, object$tags$alpha) )
  else
    cat( sprintf("%s vs %s. Test accepted: p-value: %.4e >= %.4f", object$bestMethod, object$worstMethod, object$tags$pvalue, object$tags$alpha) )
}

.testPaired <- function(bestMethod, worstMethod, statistic, pvalue, alpha, target, objetive, tags) {
  
  if (pvalue < alpha)
    outcome <- "Rejected"
  else
    outcome <- "Not Rejected"
  title <- sprintf("Wilcoxon Signed Rank Test for output \"%s\"", target)
  
  newTags <- .metaTags(title      = title,
                       target     = target,
                       objetive   = objetive,
                       alpha      = alpha,
                       statistic  = statistic,
                       pvalue     = pvalue,
                       outcome    = outcome,
                       alias      = "TestWilcoxon")
  tags <- .updateTags(tags, newTags)
  
  w <- list(
    "bestMethod"     = bestMethod,
    "worstMethod"    = worstMethod,
    "tags"           = tags
  )
  class(w) <- c("testPaired", "reportable")
  w
}
