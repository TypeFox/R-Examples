## last modified May 2008

fitted.mix <- function(object, digits = NULL, ...) 
{
    mixobj<-object
    pmat <- grpintprob(mixobj$mixdata, mixobj$parameters, mixobj$distribution, 
        mixobj$constraint)
    n <- sum(mixobj$mixdata[, 2])
    joint <- n * sweep(pmat, 2, mixobj$parameters[, 1], "*")
    mixed <- apply(joint, 1, sum)
    conditprob <- sweep(joint, 1, mixed, "/")
    if (mixobj$usecondit) {
        conditional <- sweep(conditprob, 1, apply(mixobj$mixdata[, 
            -(1:2)], 1, sum), "*")
        outlist <- list(mixed = mixed, joint = joint, conditional = conditional, 
            conditprob = conditprob)
    }
    else outlist <- list(mixed = mixed, joint = joint, conditprob = conditprob)
    if (is.null(digits)) 
        outlist
    else sapply(outlist, round, digits)
}
