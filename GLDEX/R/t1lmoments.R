`t1lmoments` <-
function (data, rmax = 4) 
{
    rmax_out <- min(rmax, 4)
    if (rmax != 4) {
        if (rmax > 4) 
            warning("The current version t1lmoments uses rmax=4.")
        rmax <- 4
    }
    data <- as.matrix(data)
    n <- dim(data)[1]
    p <- dim(data)[2]
    x <- array(, c(p, n))
    L <- array(0, c(p, rmax))
    for (j in 1:p) {
        x[j, ] <- sort(data[, j])
    }
    i <- 3:n
    s11 <- c(1, (i - 1)/(i - 2) * (n - i)/(n - i + 1))
    s12 <- c(1, (i - 1)/(i - 2) * (n - i - 1)/(n - i + 1))
    s13 <- c(1, (i - 1)/(i - 2) * (n - i - 2)/(n - i + 1))
    s14 <- c(1, (i - 1)/(i - 2) * (n - i - 3)/(n - i + 1))
    s21 <- c(1, (i - 1)/(i - 3) * (n - i)/(n - i + 1))
    s22 <- c(1, (i - 1)/(i - 3) * (n - i - 1)/(n - i + 1))
    s23 <- c(1, (i - 1)/(i - 3) * (n - i - 2)/(n - i + 1))
    s31 <- c(1, (i - 1)/(i - 4) * (n - i)/(n - i + 1))
    s32 <- c(1, (i - 1)/(i - 4) * (n - i - 1)/(n - i + 1))
    s41 <- c(1, (i - 1)/(i - 5) * (n - i)/(n - i + 1))
    s21[2] <- 1
    s22[2] <- 1
    s23[2] <- 1
    s31[2] <- 1
    s31[3] <- 1
    s32[2] <- 1
    s32[3] <- 1
    s41[2] <- 1
    s41[3] <- 1
    s41[4] <- 1
    c11 <- choose(n - 2, 1)/choose(n, 3) * cumprod(s11)
    c12 <- choose(n - 2, 2)/choose(n, 4) * cumprod(s12)
    c13 <- choose(n - 2, 3)/choose(n, 5) * cumprod(s13)
    c14 <- choose(n - 2, 4)/choose(n, 6) * cumprod(s14)
    c21 <- choose(n - 3, 1)/choose(n, 4) * cumprod(s21)
    c22 <- choose(n - 3, 2)/choose(n, 5) * cumprod(s22)
    c23 <- choose(n - 3, 3)/choose(n, 6) * cumprod(s23)
    c31 <- choose(n - 4, 1)/choose(n, 5) * cumprod(s31)
    c32 <- choose(n - 4, 2)/choose(n, 6) * cumprod(s32)
    c41 <- choose(n - 5, 1)/choose(n, 6) * cumprod(s41)
    c21[1] <- 0
    c22[1] <- 0
    c23[1] <- 0
    c31[1] <- 0
    c31[2] <- 0
    c32[1] <- 0
    c32[2] <- 0
    c41[1] <- 0
    c41[2] <- 0
    c41[3] <- 0
    for (j in 1:p) {
        L[j, 1] <- sum(c11[1:(n - 2)] * x[j, 2:(n - 1)])
        L[j, 2] <- sum((c21[1:(n - 2)] - c12[1:(n - 2)]) * x[j, 
            2:(n - 1)])/2
        L[j, 3] <- sum((c31[1:(n - 2)] - 2 * c22[1:(n - 2)] + 
            c13[1:(n - 2)]) * x[j, 2:(n - 1)])/3
        L[j, 4] <- sum((c41[1:(n - 2)] - 3 * c32[1:(n - 2)] + 
            3 * c23[1:(n - 2)] - c14[1:(n - 2)]) * x[j, 2:(n - 
            1)])/4
    }
    return(L[, 1:rmax_out])
}

