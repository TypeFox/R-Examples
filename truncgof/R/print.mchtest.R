"print.mchtest" <-
function(x, ...){
    class(x) <- "htest"
    print(x)
    cat("treshold =", x$treshold)
    cat(", simulations:", x$sim.no)
    cat("\n\n")
    invisible(x)
}
