print.ictest<-function(x,digits=4,...){
    # copy first part of  print.htest
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep = "\n")
    cat("\n")
    cat("data: ", x$data.name, "\n")
    out <- character()
    if (!is.null(x$statistic)) 
        out <- c(out, paste(names(x$statistic), "=", format(round(x$statistic, 
            4))))
     if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = digits)
        out <- c(out, paste("p-value", if (substr(fp, 1, 1) == 
            "<") fp else paste("=", fp)))
    }
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    if (!is.null(x$alt.phrase)) {
        cat("alternative hypothesis: ")
    cat(strwrap(x$alt.phrase), sep = "\n")
    }
    cat("\n")
    NU<-matrix(c(x$N,x$U),ncol=2,dimnames=list(names(x$N),
        c("n","Score Statistic*")))
    print(NU)
    if (length(x$U)>1){ cat("* like Obs-Exp, positive implies earlier failures than expected")
    } else if (x$U>0){
     cat("* postive so larger covariate values give earlier failures than expected") 
    } else if (x$U<0){
     cat("* negative so larger covariate values give later failures than expected") 
    } else cat("* zero so no trend between covariate and failure times")
    cat("\n")
 
    if (x$algorithm=="exact.mc"){
        cat(paste("p-value estimated from",x$nmc,"Monte Carlo replications"))
        cat("\n")
        cat(format(100 * attr(x$p.conf.int, "conf.level")), 
        "percent confidence interval on p-value:\n", 
            format(c(x$p.conf.int[1], x$p.conf.int[2])), "\n")
    } else if (x$algorithm=="wsr.HLY" | x$algorithm=="wsr.pclt" | x$algorithm=="wsr.mc") {
        cat(paste("p-value estimated from",x$nmc,"Monte Carlo replications"))
        cat("\n")
        if (x$algorithm=="wsr.mc"){
        cat(paste("and",x$np,"permutation resamples"))
        cat("\n")
        }
    }  

    class(x)<-"ictest"
    invisible(x)
}

