## .................................................................
is.nFactors <-
function(object) {
 if (any(class(object) == "nFactors")) return(TRUE) else return(FALSE)
 }
## .................................................................


## .................................................................
print.nFactors <-
function(x, ...) {
 if (!is.nFactors(x)) stop("Not a nFactors object")
 res <- x$nFactors
 print(res, ...)
 }
## .................................................................

## .................................................................
summary.nFactors <-
function(object, ...) {
 if (!is.nFactors(object)) stop("Not a nFactors object")
 cat("Report For a nFactors Class \n\n")
 NextMethod()
 cat(paste("Details:","\n\n"))
 print(object$detail, ...)
 cat(paste("\n\n Number of factors retained by index","\n\n"))
 print(object$nFactors)
 }
## .................................................................

