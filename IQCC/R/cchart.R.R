cchart.R <- function(x, n, type = "norm", y = NULL)
{
    if(type == "norm")
    {
        qcc(x, type = "R", xlab = "")
        resu <- alpha.risk(n)
        result <- signif(resu, 3)
        mtext(paste("Warning: Prob. of false alarm alpha = ", result," is inflated (>> 0.0027) since the normal approximation for R is not appropriated; in order to have alpha = 0.0027 the exact distribution for R must be used."), side = 1, font = 2)
    }
    if(type == "tukey")
    {
        qcc(x, type = "R", limits = c(qtukey(0.00135, n, Inf) * sd.R(y), qtukey(0.99865, n, Inf) * sd.R(y)))
    }
}