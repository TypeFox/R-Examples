gx.t.test <-
function (xbar1,s1,n1,xbar2,s2,n2) 
{
    xdif <- xbar1 - xbar2
    df <- n1 + n2 - 2
    temp1 <- ((n1-1)*s1*s1 + (n2-1)*s2*s2) / df
    temp2 <- 1/n1 + 1/n2
    t4dave <- xdif / sqrt(temp1 * temp2)
    tprob <- pt(abs(t4dave), df)
    if(tprob >= 0.9999) tprob <- 0.9999
    cat("  [xmean1, sd1, n1]:", xbar1, s1, n1, "\t[xmean2, sd2, n2]:", xbar2, s2, n2)
    cat("\n  Difference =", xdif, "\tt =", signif(t4dave,5), "\tdf =", df, 
        "\tp =", signif(tprob,4))
    if(tprob <= 0.95) { 
        cat("\n  Accept hypothesis that means are equal at the 95% confidence level\n\n") 
        }
    else {
        cat("\n  Reject hypothesis that means are equal at the 95% confidence level\n\n")
        }
    invisible()
}
