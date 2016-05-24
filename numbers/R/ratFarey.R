##
##  r a t F a r e y . R  Farey Approximation
##


ratFarey <- function(x, n, upper = TRUE) {
    sgn <- sign(x); x <- abs(x)
    xn <- trunc(x); x <- x - xn

    F1 <- c(0, 1); f1 <- F1[1]/F1[2]    # Farey series
    F2 <- c(1, 1); f2 <- F2[1]/F2[2]    #  and mediant
    repeat {
        M <- c(F1[1] + F2[1], F1[2] + F2[2]) 
        m <- M[1]/M[2] 
        if (M[2] > n) {
            if (upper) {
                F <- F2; break
            } else {
                if (x < (f1+f2)/2) {F <- F1; break}
                else               {F <- F2; break}
            }
        } else {
            if (x < m) {
                F2 <- M
                f2 <- F2[1]/F2[2]
            } else if (x > m) {
                F1 <- M
                f1 <- F1[1]/F1[2]
            } else {
                {F <- M; break}
            }
        }
    }
    F[1] <- sgn * (xn * F[2] + F[1])

    return(F)
}
