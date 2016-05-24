cchart.T2.1 <- function(T2, m, n, p)
{
    if(n < 1)
    {
        sprintf("n must be equal to or higher than 1.")
        return()
    }
    if(n == 1)
    {
        Q <- (2 * (m - 1) ^ 2) / (3 * m - 4)
        UCL <- (((m - 1) ^ 2) / m) * qbeta(1 - 0.0027, p / 2, (Q - p - 1) / 2) 
        plot(1:m, T2, main = "Hotelling T2: Individual Observations - Phase I", ylim = c(0, UCL + 1), ylab = "T2")
    }
    if(n > 1)
    {
        UCL <- ((p * (m - 1) * (n - 1)) / (m * n - m - p + 1)) * qf(1 - 0.0027, p, m * n - m - p + 1) 
        plot(1:m, T2, main = "Hotelling T2: Subgroup Observations - Phase I", ylim = c(0, UCL + 1), ylab = "T2")
    }
    lines(1:m, T2)
    mtext("UCL", side = 4, outer = F, at = UCL , padj = 0, col = 'red', font = 2)
    abline(h = UCL, lty = 2, col = 'red')
}