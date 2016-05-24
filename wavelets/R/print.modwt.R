print.modwt <- function(x, ...)
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

    cat("MODWT Wavelet Coefficients:\n")
    for(i in 1:x@level) {
        cat(paste("Level",i))
        cat("\n")
        for(j in 1:dim(x@series)[2]) {
            cat(paste(paste(paste("Series",j),":", sep = ""),""))
            printlevel(x@W[[i]][,j])
        }        
    }
    cat("\n")
  
    cat("MODWT Scaling Coefficients:\n")
    for(i in 1:x@level) {
        cat(paste("Level",i))
        cat("\n")
        for(j in 1:dim(x@series)[2]) {
            cat(paste(paste(paste("Series",j),":", sep = ""),""))
            printlevel(x@V[[i]][,j])
        }        
    }
    cat("\n")

    cat("Length of Original Series: ")
    cat(dim(x@series)[1])
    cat("\n")

    cat("Wavelet Coefficients Aligned? ")
    if(x@aligned == FALSE) {
        cat("FALSE\n")
    }
    else {
        cat("TRUE\n")
    }

    #print whether center of energy method used
    cat("Center of Energy Method Used? ")
    if(x@coe == FALSE) {
        cat("FALSE\n")
    }
    else {
        cat("TRUE\n")
    }

    cat("\n")

    cat("Boundary Method: ")
    cat(paste(toupper(substr(x@boundary, start = 1, stop =1)), substr(x@boundary, start = 2, stop = nchar(x@boundary)), sep = ""))
    cat("\n")

    cat("Number of Boundaries Coefficients per Level:\n")
    for(i in 1:length(x@n.boundary)) {
        cat(paste(paste(paste("Level",i),":", sep = ""),""))
        cat(x@n.boundary[i])
        cat("\n")  
    }
    cat("\n")  
 
    cat("Filter Class: ")
    cat(x@filter@wt.class)
    cat("\n")
    cat("Filter Name: ")
    cat(toupper(x@filter@wt.name))
    cat("\n")
    cat("Filter Length: ")
    cat(x@filter@L)
    cat("\n")
   
    invisible(x)
}
