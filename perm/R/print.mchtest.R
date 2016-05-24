print.mchtest<-function(x,...){
    class(x)<-"htest"
    print(x)
    cat(paste("p-value estimated from",x$nmc,"Monte Carlo replications"))
    cat("\n")
    cat(format(100 * attr(x$p.conf.int, "conf.level")), 
        "percent confidence interval on p-value:\n", 
            format(c(x$p.conf.int[1], x$p.conf.int[2])), "\n")

    class(x)<-"mchtest"
    invisible(x)
}

