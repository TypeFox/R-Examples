RE.ADT <- function(x){
alpha <- 0.05
n <- length(x)
s <- sd(x)
x <- sort(x)

x <- matrix(x,nrow = n)
fx <- matrix(x,nrow = n)

fx <- pnorm((x - mean(x)) / s)

i <- 1:n

S <- sum((((2 * i) - 1) / n) * (log(fx) + log(1 - fx[n + 1 - i])))
   AD2 <- - n - S
   AD2a <- AD2 * (1 + 0.75 / n + 2.25 / n ^ 2)

ifelse((AD2a >= 0.00 & AD2a < 0.200),(p <- 1 - exp( - 13.436 + 101.14 * AD2a - 223.73 * AD2a ^ 2)),
(ifelse((AD2a >= 0.200 & AD2a < 0.340),(p <- 1 - exp( - 8.318 + 42.796 * AD2a - 59.938 * AD2a ^ 2)),
(ifelse((AD2a >= 0.340 & AD2a < 0.600),(p <- exp(0.9177 - 4.279 * AD2a - 1.38 * AD2a ^ 2)),
((p <- exp(1.2937 - 5.709 * AD2a + 0.0186 * AD2a ^ 2))))))))

outList = list("Anderson-Darling Test", p = p)
    invisible(outList)

}

