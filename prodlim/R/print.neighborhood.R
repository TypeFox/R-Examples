##' @export 
"print.neighborhood" <- function(x,...){
    n <- x$n
    size <- x$size.nbh
    bw <- lapply(x$bandwidth,function(bw)round(bw,3))
    cat("Nearest neighborhoods for kernel smoothing\n\n")
    print(c(bandwidth=as.numeric(bw),kernel=x$kernel,n.obs=x$n,n.values=x$nu),quote=FALSE)
    cat("\n")
    print(c("Number of nbh's" = length(size),
            "Average size"=round(mean(size)),
            "Min size"=round(min(size)),
            "Max size"=round(max(size))))
    #  if (print.it) print(data.frame(Nbh=x$values,First=x$first.nbh,Size=size))
    invisible(x)
}
