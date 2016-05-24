## last modified May 2008

summary.mix <- function(object, digits = 4, ...) 
{
    mixobj<-object
    cat("\nParameters:\n")
    print(format(mixobj$parameters, digits = digits))
    cat("\nStandard Errors:\n")
    print(format(mixobj$se, digits = digits))
    cat("\n")
    anova.mix(mixobj)
    cat("\n")
    invisible(mixobj)
}
