cond.mode=function (Fc, method = "monoH.FC",draw=TRUE) 
{
    x = Fc$y0
    Fc = Fc$Fc
    ndist = length(x)
    if (draw) par(mfrow = c(2, 1))
    if (method == "diff") {
		if (draw) {
        plot(x, Fc, type = "l", main = "Conditional Distribution Function", 
            ylab = "Fc")}
        if (is.vector(Fc)) 
            Fc = matrix(Fc, ncol = 1)
        ind3 = 1
        der = (Fc[2:ndist, ind3] - Fc[1:(ndist - 1), ind3])/(x[2:ndist] - 
            x[1:(ndist - 1)])
        x<-x[2:(ndist)]    
        names(der) = round(x, 2)
    }
    else {
        f <- splinefun(x, Fc, method = method)
        if (draw) {curve(f(x), x[1], x[length(Fc)], main = "Conditional Distribution Function", 
            ylab = "Fc")}
        der = f(x, deriv = 1)
    }
    max.der = max(der)
    ind.mode.cond = which.max(der)
    mode.cond = x[ind.mode.cond]
#    print((x));    print((der))
	if (draw) {plot(x,der, type = "l", main = c("Conditional mode=", round(mode.cond, 
        3)), xlab = "y",ylab = "Derivative")
    abline(v = x[ind.mode.cond], col = 2, lty = 2)}
    return(list(mode.cond = mode.cond, x = x,f=der))
}

