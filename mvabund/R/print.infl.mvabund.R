print.infl.mvabund <- 
function (x, digits = max(3, getOption("digits") - 4), ...) {

    cat("Influence measures of\n\t", deparse(x$call), ":\n")
    varnames <- names(x[[1]])
    for(j in 1:length(x[[1]])){
    	cat("\n", varnames[j],":\n", sep="")

    	is.star <- apply(x$is.inf[[j]], 1, any, na.rm = TRUE)
    	print(data.frame(x$infmat[[j]], inf = ifelse(is.star, "*", " ")), 
        digits = digits, ...)
    }
    invisible(x)

}

