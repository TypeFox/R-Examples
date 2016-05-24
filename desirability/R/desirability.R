"dMax"    <- function(low, ...) UseMethod("dMax")
"dMin"    <- function(low, ...) UseMethod("dMin")
"dTarget" <- function(low, ...) UseMethod("dTarget")
"dArb"    <- function(x, ...) UseMethod("dArb")
"dBox"    <- function(low, ...) UseMethod("dBox")
"dOverall"    <- function(...) UseMethod("dOverall")
"dCategorical"    <- function(...) UseMethod("dCategorical")


dMax.default <- function(low, high, scale = 1, tol = NULL, ...)
{
   if(low >= high) stop("the low value must be greater than the high value")
   if(scale <= 0) stop("the scale parameter must be greater than zero")
      
   tmp <- list(low = low, high = high, scale = scale, missing = NA, ...)
   testSeq <- seq(tmp$low, tmp$high, length = 100)
   nonInformValue <- mean(predict.dMax(tmp, testSeq))

   structure(
      list(low = low, high = high, scale = scale, missing = nonInformValue, tol = tol, call = match.call(expand.dots = TRUE)),
      class = "dMax")      
}

dMin.default <- function(low, high, scale = 1, tol = NULL, ...)
{
   if(low >= high) stop("the low value must be greater than the high value")
   if(scale <= 0) stop("the scale parameter must be greater than zero")

   tmp <- list(low = low, high = high, scale = scale, missing = NA)
   testSeq <- seq(tmp$low, tmp$high, length = 100)
   nonInformValue <- mean(predict.dMin(tmp, testSeq))
   
   structure(
      list(low = low, high = high, scale = scale, missing = nonInformValue, tol = tol, call = match.call(expand.dots = TRUE)),
      class = "dMin")      
   
}

dTarget.default <- function(low, target, high, lowScale = 1, highScale = 1, tol = NULL, ...)
{
   if(low >= high) stop("the low value must be greater than the high value")
   if(low >= target) stop("the low value must be greater than the target")
   if(target >= high) stop("the target value must be greater than the high value")
   if(lowScale <= 0 | highScale <= 0) stop("the scale parameter must be greater than zero")

   tmp <- list(low = low, target = target, high = high, lowScale = lowScale, highScale = highScale, missing = NA)
   testSeq <- seq(tmp$low, tmp$high, length = 100)
   nonInformValue <- mean(predict.dTarget(tmp, testSeq))

   structure(
      list(low = low, target = target, high = high, lowScale = lowScale, highScale = highScale, missing = nonInformValue, tol = tol, call = match.call(expand.dots = TRUE)),
      class = "dTarget")   
}

dArb.default <- function(x, d, tol = NULL, ...)
{
   if(any(d > 1)| any(d < 0)) stop("the desirability values must be 0 <= d <= 1")
   if(length(x) != length(d)) stop("x and d must have the same length")
   if(length(x) < 2 | length(d) < 2) stop("x and d must have at least two values")
   
   ord <- order(x)
   x <- x[ord]
   d <- d[ord]
   
   tmp <- list(x = x, d = d, missing = NA)
   testSeq <- seq(min(x), max(x), length = 100)
   nonInformValue <- mean(predict.dArb(tmp, testSeq), na.rm = TRUE)

   structure(
      list(x = x, d = d, missing = nonInformValue, tol = tol, call = match.call(expand.dots = TRUE)),
      class = "dArb")   
   
}

dBox.default <- function(low, high, tol = NULL, ...)
{
   if(low >= high) stop("the low value must be greater than the high value")
      
   tmp <- list(low = low, high = high, missing = NA)
   testSeq <- seq(tmp$low, tmp$high, length = 100)
   nonInformValue <- mean(predict.dBox(tmp, testSeq))

   structure(
      list(low = low, high = high, missing = nonInformValue, tol = tol, call = match.call(expand.dots = TRUE)),
      class = "dBox")      
}

dOverall.default <- function(...)
{
   dObjs <- list(...)
   dClasses <- unlist(lapply(dObjs, class))
   if(!all(dClasses %in% c("dMax", "dMin", "dTarget", "dArb", "dBox", "dCategorical")))
      stop("some classes do not have classes in dMax, dMin, dTarget, dArb, dCategorical or dBox")
   structure(
      list(d = dObjs, call = match.call(expand.dots = TRUE)),
      class = "dOverall")   
}

dCategorical.default <- function (values, tol = NULL, ...) 
{
    if(length(values) < 2) stop("'values' should have at least two values")
    vals <- names(values)
    if(any(vals == "") | is.null(vals)) stop("'values' should be a named vector")
    if(!is.vector(values)) stop("'values' should be a named vector")
    tmp <- list(values = values, tol = tol, missing = NA, ...)
    nonInformValue <- mean(predict.dCategorical(tmp, names(values)))
    structure(list(values = values, tol = tol, missing = nonInformValue, 
        tol = tol, call = match.call(expand.dots = TRUE)), class = "dCategorical")
}
