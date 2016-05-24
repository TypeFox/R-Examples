summary.msfit <- function(object,complete=FALSE,variance=FALSE,...)
{
    if (!inherits(object, "msfit"))
        stop("'object' must be a 'msfit' object")
    cat("An object of class 'msfit'\n")
    Haz <- object$Haz
    K <- max(Haz$trans)
    msft <- unique(Haz$time) # the time points
    nt <- length(msft)
    if (nt<=12 | complete) {
        for (k in 1:K) {
            Hazk <- Haz[Haz$trans==k,]
            cat("\nTransition",k,":\n\n")
            print(Hazk,...)
        }
    }
    else {
        for (k in 1:K) {
            Hazk <- Haz[Haz$trans==k,]
            cat("\nTransition",k,"(head and tail):\n\n")
            print(head(Hazk),...)
            cat("\n...\n\n")
            print(tail(Hazk),...)
        }
    }
    if (variance) {
        cat("\n\nVariances and covariances:\n")
        varHaz <- object$varHaz
        if (nt<=12 | complete) {
            for (k1 in 1:K) {
                for (k2 in k1:K) {
                    varHazk1k2 <- varHaz[varHaz$trans1==k1 & varHaz$trans2==k2,]
                    if (k1==k2) cat("\nVariance for transition",k1,":\n\n")
                    else cat("\nCovariance between transition2",k1,"and",k2,":\n\n")
                    print(varHazk1k2,...)
                }
            }
        }
        else {
            for (k1 in 1:K) {
                for (k2 in k1:K) {
                    varHazk1k2 <- varHaz[varHaz$trans1==k1 & varHaz$trans2==k2,]
                    if (k1==k2) cat("\nVariance for transition",k1,":\n\n")
                    else cat("\nCovariance between transitions",k1,"and",k2,":\n\n")
                    print(head(varHazk1k2),...)
                    cat("\n...\n\n")
                    print(tail(varHazk1k2),...)
                }
            }
        }
    }
    return(invisible())
}
