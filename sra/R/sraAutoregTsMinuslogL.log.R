sraAutoregTsMinuslogL.log <-
function (data.mean, data.var, data.N, theor.mean, theor.var, 
    logvarME = log(1e-20)) 
{
    minuslogL <- 0
    varME <- exp(logvarME)
    if ((length(data.mean) != length(data.var)) | (length(data.mean) != 
        length(data.N)) | (length(data.mean) != length(theor.mean)) | 
        (length(data.mean) != length(theor.var))) {
        stop("vectors have different lengths")
    }
    for (t in 1:(length(theor.mean))) {
        if (is.infinite(theor.var[t]) || is.nan(theor.var[t]) || 
            (theor.var[t] < 0)) {
            return(Inf)
        }
        K <- -0.5 * data.N[t] * log(2 * pi) + 0.5 * (1 - data.N[t]) * 
            log(theor.var[t]) - 0.5 * log(theor.var[t] + data.N[t] * 
            varME)
        minuslogL <- minuslogL - K + data.N[t] * data.mean[t] + 
            data.N[t] * (theor.var[t] * (data.var[t] + (data.mean[t] - 
                theor.mean[t])^2) + data.N[t] * data.var[t] * 
                varME)/(2 * theor.var[t] * (theor.var[t] + data.N[t] * 
                varME))
    }
    return(minuslogL)
}
