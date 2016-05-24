#' @S3method print treemox
print.treemox <-
function(x, ...)
{
    cat("SEGMENTATION TREES IN PLS PATH MODELING", "\n")
    cat("----------------------------------------------------------", "\n")    
    cat("Tree Specification:", "\n")
    cat("1  Mox Algorithm       ", x$model$mox, "\n")
    cat("2  Threshold signif    ", x$model$signif, "\n")
    if (x$model$size<1) {
        cat("3  Node size limit(%)  ",  x$model$size, "\n")
    } else {
        cat("3  Node size limit(#)  ",  x$model$size, "\n")
    }
    cat("4  Tree depth level    ",  x$model$deep, "\n")
    cat("\n")
    cat("----------------------------------------------------------", "\n")    
    cat("Segmentation Variables:", "\n")
    print(x$model$df.exev, print.gap=2)
    cat("----------------------------------------------------------", "\n\n")    
    cat("$MOX", "\n")
    print(x$MOX, digits=3, print.gap=2)
    invisible(x)
}

