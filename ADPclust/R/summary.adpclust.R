##' Summarizes results from adpclust function from the ADPclust package.
##'
##' @title Summary of adpclust
##' @param object object with class "adpclust", return of adpclust()
##' @param ... other arguments. NOT used in theis version. 
##' @return NULL
##' @author Yifan Xu
##' @export

summary.adpclust <- function(object, ...){
    cat("-- ADPclust Procedure -- \n\n")
    cat("Number of variables: \t", attr(object$dat, "p"), "\n")
    cat("Number of obs.: \t", attr(object$dat, "n"), "\n")
    cat("Centroids selection: \t", ifelse(object$testPars$centroids == "auto", "Automatic", "Manual"), "\n")
    cat("Bandwith selection: \t", paste0(object$testPars$htype, " (", round(object$h,2), ")"),"\n")    
    cat("Number of clusters: \t", length(object$centers), "\n")
    cat("Avg. Silhouette: \t", object$score, "\n")
    cat("\nf(x): \n")
    print(summary(object$f))
    cat("\ndelta(x): \n")
    print(summary(object$delta))
    invisible(NULL)
}

