rMVN <- function(n, mean, sigma) {
            L  <- mroot(sigma)
            lL <- ncol(L)
            t(mean + L%*%matrix(rnorm(lL*n), lL, n)) 
}

