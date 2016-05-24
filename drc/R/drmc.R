"drmc" <- function(constr = FALSE, errorm = TRUE, maxIt = 500, method = "BFGS", 
noMessage = FALSE, relTol = 1e-7, rmNA = FALSE, useD = FALSE, trace = FALSE, 
otrace = FALSE, warnVal = -1, dscaleThres = 1e-15, rscaleThres = 1e-15)
{
    return(list(
                constr = constr,
                errorm = errorm,
                maxIt = maxIt, 
                method = method,
                noMessage = noMessage,
                relTol = relTol,
                rmNA = rmNA, 
                useD = useD,
                trace = trace,
                otrace = otrace,
                warnVal = warnVal,
                dscaleThres = dscaleThres,
                rscaleThres = rscaleThres
                ))
}
