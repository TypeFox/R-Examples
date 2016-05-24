#Internal Function

zeta.fun <- function (x) 
{
    if (any(x < 1)) 
        stop("Invalid input for Riemann's zeta function.", "\n")
    a = 12
    k = 8
    B = c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510)
    ans = 0
    ans = sum(1/c(1:(a - 1))^x)
    ans = ans + 1/((x - 1) * a^(x - 1)) + 1/(2 * a^x)
    term = (x/2)/a^(x + 1)
    ans = ans + term * B[1]
    for (mm in 2:k) {
        term = term * (x + 2 * mm - 2) * (x + 2 * mm - 3)/(a * 
            a * 2 * mm * (2 * mm - 1))
        ans = ans + term * B[mm]
    }
    ans
}