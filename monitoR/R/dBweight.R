# From seewave?

dBweight <-
function (f, dBref = NULL) 
{
    num <- (12200^2 * f^4)
    den <- (f^2 + 20.6^2) * sqrt((f^2 + 107.7^2) * (f^2 + 737.9^2)) * 
        (f^2 + 12200^2)
    A <- 2 + 20 * log10(num/den)
    num <- (12200^2 * f^3)
    den <- (f^2 + 20.6^2) * sqrt((f^2 + 158.5^2)) * (f^2 + 12200^2)
    B <- 0.17 + 20 * log10(num/den)
    num <- (12200^2 * f^2)
    den <- (f^2 + 20.6^2) * (f^2 + 12200^2)
    C <- 0.06 + 20 * log10(num/den)
    a <- f/(6.8966888496476 * 10^-5)
    h <- ((1037918.48 - f^2)^2 + (1080768.16 * f^2))/((9837328 - 
        f^2)^2 + (11723776 * f^2))
    b <- sqrt(h/((f^2 + 79919.29) * (f^2 + 1345600)))
    D <- 20 * log10(a * b)
    result <- list(A = A, B = B, C = C, D = D)
    if (!is.null(dBref)) 
        result <- list(A = dBref + A, B = dBref + B, C = dBref + 
            C, D = dBref + D)
    return(result)
}
