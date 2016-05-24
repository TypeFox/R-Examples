##
## Unit testing of Linear Program with cone constraints
test.LPCC <- function(){
    ## Creating LP
    q <- c(-6, -4, -5)
    ## NNO constraint
    G <- matrix(c(16, -14, 5,
                  7, 2, 0),
                nrow = 2, ncol = 3, byrow = TRUE)
    h <- c(-3, 5)
    nno1 <- nnoc(G = G, h = h)
    ## SOC constraint
    F1 <- matrix(c(8, 13, -12,
                   -8, 18, 6,
                   1, -3, -17),
                 nrow = 3, ncol = 3, byrow = TRUE)
    g1 <- c(-2, -14, -13)
    d1 <- c(-24, -7, 15)
    f1 <- 12
    soc1 <- socc(F = F1, g = g1, d = d1, f = f1)
    ## PSD constraint
    F1 <- matrix(c(7, -5, 1,
                   -5, 1, -7,
                   1, -7, -4),
                 nrow = 3, ncol = 3, byrow = TRUE)
    F2 <- matrix(c(3, 13, -6,
                   13, 12, -10,
                   -6, -10, -28),
                 nrow = 3, ncol = 3, byrow = TRUE)
    F3 <- matrix(c(9, 6, -6,
                   6, -7, -7,
                   -6, -7, -11),
                 nrow = 3, ncol = 3, byrow = TRUE)
    F0 <- matrix(c(68, -30, -19,
                   -30, 99, 23,
                   -19, 23, 10),
                 nrow = 3, ncol = 3, byrow = TRUE)
    psd1 <- psdc(Flist = list(F1, F2, F3), F0 = F0)
    ## Using main function of package
    ans <- cccp(q = q, cList = list(nno1, soc1, psd1))
    checkTrue(ans$status == "optimal")
    return()
}
