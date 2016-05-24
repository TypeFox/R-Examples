na.trim <- function(object, ...) UseMethod("na.trim")
na.trim.default <- function (object, sides = c("both", "left", "right"), 
	is.na = c("any", "all"), ...)
{
   is.na <- match.arg(is.na)
   nisna <- if (is.na == "any" || length(dim(object)) < 2)  {
	complete.cases(object)
   } else rowSums(!is.na(object)) > 0
   idx <- switch(match.arg(sides), left = cumsum(nisna) > 0,
       right = rev(cumsum(rev(nisna) > 0)>0),
       both = (cumsum(nisna) > 0) & rev(cumsum(rev(nisna)) > 0))
   if (length(dim(object)) < 2)
       object[idx]
   else
       object[idx,, drop = FALSE]
}

## need a 'ts' method because indexing destroys ts attributes
na.trim.ts <- function (object, ...)
{
    as.ts(na.trim(as.zoo(object), ...))
}
