##
##  p o l y l o g . R
##


polylog <- function(z, n) {
    stopifnot(is.numeric(z), is.numeric(n))
    if (length(n) != 1 || floor(n) != ceiling(n) || n < -4)
        stop("Argument 'n' must be a natural number n >= -4 .")
    if (length(z) > 1 || abs(z) >= 1)
        stop("Argument 'z' must be a single real number with abs(z) < 1 .")

    b <- function(i) zeta_(n - i)
    S <- function(n, z, j) {
        out <- 0
        for (k in 1:j) out <- out + z^k/k^n
        return(out)
    }
    eta_ <- function(x, j) {
        out <- 0
        for (k in 1:j) out <- out + (-1)^(k+1) / k^x
        return(out)
    }
    zeta_ <- function(x) {
        prefactor <- 2^(x-1) / ( 2^(x-1)-1 )
        numerator <- 1 + 36*2^x*eta_(x,2) + 315*3^x*eta_(x,3) + 1120*4^x*eta_(x,4) + 
                      1890*5^x*eta_(x,5) + 1512*6^x*eta_(x,6) + 462*7^x*eta_(x,7)
        denominator <- 1 + 36*2^x + 315*3^x + 1120*4^x + 1890*5^x + 1512*6^x + 
                      462*7^x
        return(prefactor * numerator / denominator)
    }

    alpha <- -log(z)
    
    if (abs(z) > 0.55) {
        preterm <- gamma(1-n)/alpha^(1-n)
        nominator <- b(0) - 
              alpha * ( b(1) - 4*b(0)*b(4)/7/b(3) ) + 
              alpha^2 * ( b(2)/2 + b(0)*b(4)/7/b(2) - 4*b(1)*b(4)/7/b(3) ) - 
              alpha^3 * ( b(3)/6 - 2*b(0)*b(4)/105/b(1) + b(1)*b(4)/7/b(2) - 
                          2*b(2)*b(4)/7/b(3) )
        denominator <- 1 + alpha*4*b(4)/7/b(3) + 
              alpha^2 * b(4) / 7 / b(2) + 
              alpha^3 * 2 * b(4) / 105 / b(1) + 
              alpha^4 * b(4) / 840 / b(0)
        y <- preterm + nominator / denominator

    } else {
        nominator <- 6435 * 9^n * S(n,z,8) - 27456 * 8^n * z*S(n,z,7) + 
              48048 * 7^n * z^2 * S(n,z,6) - 44352 * 6^n * z^3 * S(n,z,5) + 
              23100 * 5^n * z^4 * S(n,z,4) -  6720 * 4^n * z^5 * S(n,z,3) + 
              1008 * 3^n * z^6 *S(n,z,2) - 64 * 2^n * z^7 * S(n,z,1)
        denominator <- 6435 * 9^n - 27456 * 8^n * z + 
              48048 * 7^n * z^2 - 44352 * 6^n * z^3 + 
              23100 * 5^n * z^4 -  6720 * 4^n * z^5 + 
              1008 * 3^n * z^6 - 64 * 2^n * z^7 + 
              z^8
        y <- nominator / denominator
    }

    return(y)
}

