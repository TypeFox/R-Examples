print.wt.filter <- function(x, ...)
{
    printlevel <- function(coef.level) {
        if(length(coef.level) > 6) {
            coef.level.begin <- formatC(coef.level[1:3], format = "e")
            coef.level.end <- formatC(coef.level[(length(coef.level)-2):(length(coef.level))], format = "e")
            cat(c(coef.level.begin , "..." , coef.level.end))
            cat("\n") 
        }
        else {
            cat(formatC(coef.level, format = "e"))
            cat("\n")   
        }
    }

    cat("Filter Class: ")
    cat(x@wt.class)
    cat("\n")
    cat("Name: ")
    cat(toupper(x@wt.name))
    cat("\n")
    cat("Length: ")
    cat(x@L)
    cat("\n")
    cat("Level: ")
    cat(x@level)
    cat("\n")
    cat("Wavelet Coefficients: ")
    printlevel(x@h)
    cat("Scaling Coefficients: ")
    printlevel(x@g)
    cat("\n")

    invisible(x)
}

