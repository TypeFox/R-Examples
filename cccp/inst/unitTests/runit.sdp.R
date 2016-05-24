##
## Unit testing of Linear Program with PSD constraints
test.SDP <- function(){
    ## Creating SDP
    ## Objective
    q <- c(1, -1, 1)
    ## First PSD cone constraint
    F1 <- matrix(c(-7, -11, -11, 3), nrow = 2, ncol = 2)
    F2 <- matrix(c(7, -18, -18, 8), nrow = 2, ncol = 2)
    F3 <- matrix(c(-2, -8, -8, 1), nrow = 2, ncol = 2)
    F0 <- matrix(c(33, -9, -9, 26), nrow = 2, ncol = 2)
    psd1 <- psdc(Flist = list(F1, F2, F3), F0 = F0)
    ## Second PSD cone constraint
    F1 <- matrix(c(-21, -11, 0, -11, 10, 8, 0, 8, 5), nrow = 3, ncol = 3)
    F2 <- matrix(c(0, 10, 16, 10, -10, -10, 16, -10, 3), nrow = 3, ncol = 3)
    F3 <- matrix(c(-5, 2, -17, 2, -6, 8, -17, 8, 6), nrow = 3, ncol = 3)
    F0 <- matrix(c(14, 9, 40, 9, 91, 10, 40, 10, 15), nrow = 3, ncol = 3)
    psd2 <- psdc(Flist = list(F1, F2, F3), F0 = F0)
    ## Using main function of package
    ans <- cccp(q = q, cList = list(psd1, psd2))
    checkTrue(ans$status == "optimal")
    return()
}
