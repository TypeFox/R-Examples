print.gnsc <-
function(x, ...){
    cat("Group Nearest Shrunken Centroids Output\n")
    cat("Number of Threshold Values:",x$nlambda, "\n")
    out = cbind(x$lambda, x$nonzero, x$errors)
    colnames(out) = c("threshold", " nonzero", "errors")
    print(out)
    #cat("Threshold Values:", x$lambda, "\n")
    #cat("Number of Nonzero Variable Groups:", x$nonzero, "\n")
    #cat("estimated errors:", x$errors,  "\n")
    cat("predicted sample classes: \n")
    print(x$yhat)
}
