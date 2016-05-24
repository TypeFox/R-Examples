`SiddiquiMatrix` <-
function(phi)
{
    p <- length(phi)
    phis <- c(-1, phi)
    A <- matrix(numeric(p^2), nrow = p, ncol = p)
    for(j in 1:p) 
        for(i in 1:p) 
            if(j > i) 
                A[i, j] <- A[j, i]
            else {
                k <- 1:min(i, j)
                A[i, j] <- sum(phis[1 + i - k] * phis[1 + j - k] - 
                    phis[1 + p + k - i] * phis[1 + p + k - j])
            }
    A
}

