"dplot" <-
function(x, distn, parm, H = NA, verticals = FALSE, ...) 
{
    if (!is.function(try(get(distn), silent = TRUE)))
       stop("'distn' must be a character of a distribution function")
    if (is.na(H)) {
       H <- -Inf
       warning("no treshold value specified")}
    dots <- list(...)
    pdist <- get(distn, mode = "function")
    qdist <- get(paste("q", substring(distn, 2), sep = ""), mode = "function")

    Fn <- edf(x, distn, parm, H)
    f  <-  function(x) do.call("pdist", c(list(x), parm))
         
    xv <- knots(Fn)
    rx <- range(xv)
    dr <- max(0.08 * diff(rx), median(diff(xv)))
    xlim <- rx + dr * c(-1, 1)
    ylim <- c(0,1)
       
    ti <- list(xlab = deparse(substitute(x)), ylab = "F, Fn",
        main = "Distribution Function", type = "l", lwd = 2, col = "slateblue")
        
    nm <- match(names(ti), names(dots))
    na <- !is.na(nm)  
    ti[na] <- dots[nm[na]] 
    dots <- dots[-nm[na]]
    
    do.call("curve", c(list(as.name("f"), xlim = xlim, ylim = ylim), ti, dots))
    par(new = TRUE)
    plot.stepfun(Fn, verticals = verticals, xlab = "", ylab = "", main = "")
    abline(h = c(0, 1), col = "gray70", lty = 2) 
    do.call("title", ti)
}
