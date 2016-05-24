"print.MCtest" <-
function(x,...){
    cat(paste("MCtest using MCbound of\n    type=",x$type, "with parms=\n"))
    cat("\n")
    print(x$parms,...)
    cat("\n")
    cat(paste("p-value=",format(x$p.value)))
    cat("\n")
    cat(format(100 * attr(x$p.value.ci, "conf.level")), "percent confidence interval (on p-value):\n", 
            format(c(x$p.value.ci[1], x$p.value.ci[2])), "\n")
    invisible(x)
}

