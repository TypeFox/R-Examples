ar1est0 <-
function (z) 
{
    n <- length(z)
    stopifnot(n > 3)
    a <- sum(z^2)
    b <- sum(z[-1] * z[-n])
    c <- sum(z[c(-1, -n)]^2)
    i <- complex(1, 0, 1)
    x <- ((-16) * b^3 + 18 * a * b * c + 24 * b^3 * n - 
            27 * a * b * c * n - 9 * b * c^2 * n - 12 * b^3 * 
            n^2 + 9 * a * b * c * n^2 + 27 * b * c^2 * n^2 + 
            2 * b^3 * n^3 - 18 * b * c^2 * n^3)
    y <- (-b^2 * (-2 + n)^2 + 3 * c * (-1 + n) * (-a - c * 
            n))
    f <- complex(1, x^2 + 4 * y^3, 0)
    z <- (x + sqrt(f))^(1/3)
    g <- x^2 + 4 * y^3
    z1 <- (x + (-g)^(1/2) * i)^(1/3)
    part1 <- (n - 2) * b/(3 * c * (n - 1))
    part2 <- (1 - sqrt(3) * i) * y/(3 * 2^(2/3) * c * (n - 
            1) * z)
    part3 <- (1 + sqrt(3) * i) * z/(6 * 2^(1/3) * c * (n - 
            1))
    Re(part1 + part2 - part3)
}
