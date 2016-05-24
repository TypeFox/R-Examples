`dlvcomp2` <-
function (N, alpha, rd = c(1, 1)) 
{
    N1.t1 <- N[1] + rd[1] * N[1] * (1 - alpha[1, 1] * N[1] - 
        alpha[1, 2] * N[2])
    N2.t1 <- N[2] + rd[2] * N[2] * (1 - alpha[2, 1] * N[1] - 
        alpha[2, 2] * N[2])
    c(N1.t1, N2.t1)
}
