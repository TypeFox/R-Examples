## .................................................................
summary.nScree <-
function(object, ...) {
 if (!is.nScree(object)) stop("Not a nScree object")
 cat("Report For a nScree Class \n\n")
 #digits <- 2
 NextMethod()
 cat(paste("Details:",object$Model,"\n\n"))
 object$Analysis[,c(1:5,7)] <- round(object$Analysis[,c(1:5,7)], ...)
 print(object[[2]])
 cat(paste("\n\n Number of factors retained by index","\n\n"))
 print(object[[1]])
 }
## .................................................................


## .................................................................
print.nScree <-
function(x, ...) {
 res <- x[[1]]
 print(res, ...)
 }
## .................................................................


## .................................................................
plot.nScree <-
function(x, ...) {
 plotnScree(x, ...)
 }
## .................................................................


## .................................................................
is.nScree <-
function(object) {
 if (class(object) == "nScree") return(TRUE) else return(FALSE)
 }
## .................................................................
