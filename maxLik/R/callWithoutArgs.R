## strip arguments "args" and call the function with name "fName" thereafter
callWithoutArgs <- function(theta, fName, args, ...) {
   f <- match.call()
   f[ args ] <- NULL
   f[[1]] <- as.name(fName)
   names(f)[2] <- ""
   f[["fName"]] <- NULL
   f[["args"]] <- NULL
   f1 <- eval(f, sys.frame(sys.parent()))
   return( f1 )
}
