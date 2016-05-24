LSP <-
function (N, LTPD, beta, p = seq(0, 0.3, 0.001), Plots = TRUE) 
{
    f = 1 - (beta^(1/(LTPD * N)))
    n = round(f * N)
    OC = (1 - f)^(N * p)
    AOQ = (N - n) * p * OC/N
    ATI = n * OC + N * (1 - OC)
    results = list(p = p, OC = OC, n = rep(n, length(p)), AOQ = AOQ, 
        ATI = ATI)
    class(results) = "AccSampPlan"
    if (Plots) {
        par(mfrow = c(2, 2))
        plot(OC ~ p, type = "l", ylab = "Probability of Acceptance", 
            xlab = "Fraction Nonconforming p")
        title(paste("f = ", formatC(f)))
        plot(rep(n, length(p)) ~ p, type = "l", ylab = "Uncurtailed sample size", 
            xlab = "Fraction Nonconforming p")
        title(paste("n = ", formatC(n)))
        plot(AOQ ~ p, type = "l", ylab = "AOQ", xlab = "Fraction Nonconforming p")
        title(paste("AOQL = ", formatC(max(AOQ))))
        plot(ATI ~ p, type = "l", ylab = "ATI", xlab = "Fraction Nonconforming p")
    }
    return(results)
}
