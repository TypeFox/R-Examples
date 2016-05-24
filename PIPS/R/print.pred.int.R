#####################################################################
# PROGRAM:       print.pred.int.R
# PURPOSE:       Function to print objects of class "pred.int".
# INPUT:         
#   x:  Object of type pred.int
#   pi.count: Number of predicted intervals to print
#   digits:   Digits used when printing confidence intervals
# OUTPUT:        pred.int object printed to usual location
# MACROS USED:   None
# CALLED BY:     None
# AUTHOR:        Daniel Muenz
# CREATION DATE: 2010
# NOTES: 
# MODIFICATIONS:
# [RG20120107] Ray Griner standardized stop message.
#####################################################################
print.pred.int <- function(x, pi.count=8,
                           digits=max(3, getOption("digits")-3), ...)
{
    ## Make sure pi.count is a number >= 2.
    if (!is.numeric(pi.count) || length(pi.count)!=1 || pi.count<2)
        stop("'pi.count' must be greater than 1") ## [RG20120107]

    orig.x <- x
    if (is.null(names(x$obs.n))) {
        names(x$obs.n) <- " "
        names(x$ci) <- " "
        names(x$pi) <- " "
    }
    comp.names <- names(x$ci)

    ## Print out the sample sizes.
    cat("Sample sizes:\n")
    print( cbind(Observed=x$obs.n,
                 Simulated=x$sim.n,
                 Total=x$obs.n+x$sim.n), ... )

    ## Print out the confidence intervals.
    cat("\nPoint estimates and ", 100*x$obs.conf.level,
        "% confidence intervals from observed data:\n", sep="")
    cis <- matrix(0, nrow=length(x$ci), ncol=3,
                  dimnames=list(comp.names,
                                c("Point","Lower Bound","Upper Bound")))
    for(comp in comp.names)
        cis[comp,] <- x$ci[[comp]]
    print( cis, digits=digits, ... )

    ## Print out the predicted intervals.
    cat("\nPoint estimates and ", 100*x$conf.level,
        "% predicted intervals from observed+simulated data:\n", sep="")
    for(comp in comp.names) {
        if (comp != " ") cat(comp, ":\n", sep="")
        pis <- x$pi[[comp]]
        dimnames(pis)[[2]] <- c("Point","Lower Bound","Upper Bound")
        ## Make sure we only show as many predicted intervals as is
        ## specified by pi.count.
        if ((iters <- nrow(pis)) > pi.count) {
            block.size <- as.integer(pi.count/2)
            pis <- pis[c(1:(block.size+1), (iters-(pi.count-block.size)+1):iters),]
            dimnames(pis)[[1]][block.size+1] <- "..."
            pis[block.size+1,] <- c(NA, NA, NA)
        }
        print( pis, na.print=" ", digits=digits, ... )
        cat("\n")
    }

    ## Return a copy of the object passed.  This seems to be a
    ## convention.
    invisible(orig.x)
}
