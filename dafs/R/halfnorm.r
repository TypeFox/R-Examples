halfnorm = function(fit){
    n = length(residuals(fit))
    p = ppoints(2*n+1)
    alpha = p[p>0.5]

    absRes = abs(residuals(fit))
    resLabs = as.character(1:n)
    o = order(absRes)
    absRes = absRes[o]
    resLabs = resLabs[o]

    plot(qnorm(alpha), absRes, xlab = 'Standard normal quantiles',
         ylab = 'abs(Res)')

    big = which(alpha>0.95)
    if(length(big)>0)
        text(qnorm(alpha[big]),absRes[big],resLabs[big],
             cex = 0.7, adj = c(1,0))
}
