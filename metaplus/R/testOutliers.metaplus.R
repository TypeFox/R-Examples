testOutliers <-
  ## Short form for generic function 
  function(object,R=999) UseMethod("testOutliers")

testOutliers.metaplus <- function(object,R=999) {
  if (!inherits(object, "metaplus"))
    stop("Use only with 'metaplus' objects.\n")

  if (object$justfit) stop("Cannot use with objects fitted with justfit=TRUE")
  
  if (object$random=="normal") stop("cant test for outliers with normal random effects model")
  testOutliers <- switch(object$random,
                "t-dist"=testOutliers.profilet.metaplus(object,R),
                "mixture"=testOutliers.profilemix.metaplus(object,R))
  class(testOutliers) <- "testOutliers"
  return(testOutliers)
}

summary.testOutliers <- function(object,...) {
  if (!inherits(object, "testOutliers"))
    stop("Use only with 'testOutliers' objects.\n")
  thesummary <- list(pvalue=object$pvalue,observed=object$observed)
  class(thesummary) <- "summary.testOutliers"
  thesummary
}

print.summary.testOutliers <- function(x,...) {
  if (!inherits(x, "summary.testOutliers"))
    stop("Use only with 'summary.testOutliers.metaplus' objects.\n")
  cat("Observed LRT statistic ")
  cat(sprintf("%1.1f",x$observed))
  cat(" p value ")
  thepvalue <- format.pval(x$pvalue, digits = 4, eps = 1.0e-4) 
  cat(sprintf("%s\n",thepvalue))
}
