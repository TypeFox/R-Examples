## last modified May 2008

print.mix <- function(x, digits = 4, ...) 
{
    mixobj<-x
    constr <- mixobj$constraint
    param <- mixobj$parameters
    cat("\nParameters:\n")
    if (constr$conpi == "PFX") 
        param <- cbind(param, fixpi = constr$fixpi)
    if (constr$conmu == "MFX") 
        param <- cbind(param, fixmu = constr$fixmu)
    if (constr$consigma == "SFX") 
        param <- cbind(param, fixsigma = constr$fixsigma)
    if (constr$consigma == "BINOM") 
        param <- cbind(param, size = constr$size)
    if (constr$consigma == "NBINOM") 
        param <- cbind(param, size = constr$size)
    print(format(param, digits = digits))
    cat("\nDistribution:\n")
    print(mixobj$distribution)
    cat("\nConstraints:\n")
    if (constr$consigma == "FCV") 
        print(c(conpi = constr$conpi, conmu = constr$conmu, consigma = constr$consigma, 
            cov = constr$cov))
    else print(c(conpi = constr$conpi, conmu = constr$conmu, 
        consigma = constr$consigma))
    cat("\n")
    invisible(mixobj)
}
