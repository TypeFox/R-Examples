"lillie.test" <-
function (x) 
{
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    if (n < 5) 
        stop("sample size must be greater than 4")
    p <- pnorm((x - mean(x))/sd(x))
    Dplus <- max(seq(1:n)/n - p)
    Dminus <- max(p - (seq(1:n) - 1)/n)
    K <- max(Dplus, Dminus)
    if (n <= 100) {
        Kd <- K
        nd <- n
    }
    else {
        Kd <- K * ((n/100)^0.49)
        nd <- 100
    }
    pvalue <- exp(-7.01256 * Kd^2 * (nd + 2.78019) + 2.99587 * 
        Kd * sqrt(nd + 2.78019) - 0.122119 + 0.974598/sqrt(nd) + 
        1.67997/nd)
    if (pvalue > 0.1) {
        KK <- (sqrt(n) - 0.01 + 0.85/sqrt(n)) * K
        if (KK <= 0.302) {
            pvalue <- 1
        }
        else if (KK <= 0.5) {
            pvalue <- 2.76773 - 19.828315 * KK + 80.709644 * 
                KK^2 - 138.55152 * KK^3 + 81.218052 * KK^4
        }
        else if (KK <= 0.9) {
            pvalue <- -4.901232 + 40.662806 * KK - 97.490286 * 
                KK^2 + 94.029866 * KK^3 - 32.355711 * KK^4
        }
        else if (KK <= 1.31) {
            pvalue <- 6.198765 - 19.558097 * KK + 23.186922 * 
                KK^2 - 12.234627 * KK^3 + 2.423045 * KK^4
        }
        else {
            pvalue <- 0
        }
    }
    RVAL <- list(statistic = c(D = K), p.value = pvalue, method = "Lilliefors (Kolmogorov-Smirnov) normality test", 
        data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}
