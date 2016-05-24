cchart.T2.2 <- function(T2II, m, n, j, t, p, datum = NULL, stats = NULL, T2 = NULL)
{
    if(n == 1)
        UCL <- ((p * (m + 1) * (m - 1)) / ((m ^ 2) - m * p)) * qf(1 - 0.0027, p, m - p) 
    if(n > 1)
        UCL <- ((p * (m + 1) * (n - 1))/(m * n - m - p + 1)) * qf(1 - 0.0027, p, m * n - m - p + 1)

    old = FALSE
    if(is.null(T2) == FALSE)
    {
        plot(c(1:m, j + m + 1), c(T2, T2II[1]), ylim = c(0, UCL + 1), xlim = c(1, t + m + 1), ylab = "T2", xlab = "Sample", pch = 16, xaxt = 'n')
        old = TRUE
    }
    else
    {
        if(is.null(T2) && is.null(stats) == FALSE)
        {
            T2 <- T2.1(stats, m, n)
            plot(c(1:m, j + m + 1), c(T2, T2II[1]), ylim = c(0, UCL + 1), xlim = c(1, t + m + 1), ylab = "T2", xlab = "Sample", pch = 16, xaxt = 'n')
            old = TRUE
	  }
        else
        {
            if(is.null(T2) && is.null(stats) && is.null(datum) == FALSE)
            {
                stats <- stats(datum, m, n, p)
                T2 <- T2.1(stats, m, n)
                plot(c(1:m, j + m + 1), c(T2, T2II[1]), ylim = c(0, UCL + 1), xlim = c(1, t + m + 1), ylab = "T2", xlab = "Sample", pch = 16, xaxt = 'n')
                old = TRUE
		}
            else
                if(is.null(T2) && is.null(stats) && is.null(datum))
                    plot(j, T2II[1], ylim = c(0, UCL + 1), xlim = c(1, t), ylab = "T2", xlab = "Sample", pch = 16)
        }
    }
    if(n == 1)
        title("Hotelling T2: Individual Observations - Phase II")
    if(n > 1)
        title("Hotelling T2: Subgroup Observations - Phase II")
    
    mtext("UCL", side = 4, outer = F, at = UCL , padj = 0, col = 'red', font = 2)
    abline(h = UCL, lty = 2, col = 'red')
    if(old == TRUE)
    {
        axis(1, at = 1:m, labels = 1:m)
        axis(1, at = (m+2):(m+t+1), labels = 1:t)
        lines(T2)
        abline(v = m + 1)
    }
}