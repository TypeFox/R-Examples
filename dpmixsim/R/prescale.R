# data scaling
prescale <-
function(xv)
{
     ## cat("data vector length: ", length(xv),"\n")
     ## xv <- sort(xv) # ! WITH SORT ?
     minxv <- min(xv)
     maxxv <- max(xv)
     x <- (xv-minxv) / (maxxv- minxv)
     invisible(x)
}


