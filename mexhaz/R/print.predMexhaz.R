print.predMexhaz <- function(x, ...){
    cat("Results:\n")
    print(head(x$results))
    cat(paste("... (dimensions: ",dim(x$results)[1]," rows and ",dim(x$results)[2]," columns)\n",sep=""))
}
