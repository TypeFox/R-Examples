summary.dwt <- function(object, ...)
{
    cat("Number of Points in Original Series: ")
    cat(dim(object@series)[1])
    cat("\n")

    cat("Number of Levels Decomposed: ")
    cat(object@level)
    cat("\n")

    cat("Filter Name: ")
    cat(toupper(object@filter@wt.name))
    cat("\n")

    cat("Boundary Method: ")
    bmethod <- object@boundary
    ncn <- nchar(bmethod)
    cat(paste(toupper(substr(object@boundary, start = 1, stop =1)), substr(object@boundary, start = 2, stop = nchar(object@boundary)), sep = ""))
    cat("\n")

    cat("Sum of Squares of Wavelet Coefficients:\n")
    for(i in 1:object@level) {
        cat(paste("Level",i))
        cat("\n")
        for(j in 1:dim(object@series)[2]) {
            cat(paste(paste(paste("Series",j),":", sep = ""),""))
            cat(sum(object@W[[i]][,j]^2))
            cat("\n")
        }        
    }
    cat("\n")

    cat("Sum of Squares of Scaling Coefficients:\n")
    for(i in 1:object@level) {
        cat(paste("Level",i))
        cat("\n")
        for(j in 1:dim(object@series)[2]) {
            cat(paste(paste(paste("Series",j),":", sep = ""),""))
            cat(sum(object@V[[i]][,j]^2))
            cat("\n")
        }        
    }
    cat("\n")

    invisible(object)
}

