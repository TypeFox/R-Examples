print.fptl <-
function (x, ...) 
{
    if (!is.fptl(x)) 
        stop(paste(sQuote("x"), "is not of class", shQuote("fptl")))

    cutoff <- options()$deparse.cutoff
    options(deparse.cutoff = 175)

    
    cat("\nAn object of class", shQuote("fptl"), "containing")
    cat("\n   $ x: a sequence of ", length(x$x), " values from ", 
        format(x$x[1], ...), " to ", format(x$x[length(x$x)], ...), sep = "")
    cat("\n   $ y: the values of the FPTL function on x")
    cat("\n\nCall:\n")
    print(attr(x, "Call"))
 
    dp <- attr(x, "Call")[["dp"]]
    if (is.name(dp)) cat("\n", dp, ":", sep="") else cat("\ndp:")
    print(attr(x, "dp"))
     
    v <- attr(x, "vars")
    if (!is.null(v)){
		cat("Values of names which occur in Call (except in dp argument):\n")
		print(v)		
    }

    options(deparse.cutoff = cutoff)
}
