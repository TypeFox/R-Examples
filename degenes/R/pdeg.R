pdeg <-
function (treatment, control, alpha = 0.05) 
{
    cat(date(), fill = TRUE)
    n <- dim(treatment)[1]
    J1 <- dim(treatment)[2]
    J2 <- dim(control)[2]
    m1 <- J1/2
    m2 <- J2/2
    X1 <- (treatment[, 1:m1] - treatment[, (m1 + 1):J1])/2
    X2 <- (control[, 1:m2] - control[, (m2 + 1):J2])/2
    Z.big <- z.b(treatment, control, X1, X2, m1, m2)
    z.small <- z.s(X1, X2, m1, m2)
    new1 <- Z.big
    new2 <- z.small
    z <- seq(1.5 * min(min(new1), min(new2)), 1.5 * max(max(new1), 
        max(new2)), 0.01)
    lz <- length(z)
    integral.factor <- 0.005
    kern1 <- kern(new1, z, lz)
    kern2 <- kern(new2, z, lz)
    if (max(kern2) <= max(kern1)) {
        print("No differentially genes.", quote = FALSE)
    }
    if (max(kern2) > max(kern1)) {
        p.value <- rep(0, n)
        c <- 0
        n2 <- length(new1)
        band1 <- (sqrt(var(new1)) * (n2^(-1/5))) * 1.144
        band2 <- (sqrt(var(new2)) * (n2^(-1/5))) * 1.144
        Tk.hat <- kern2[which(kern1 > 0)]/kern1[which(kern1 > 
            0)]
        for (i in 1:n) {
            zkk <- rep(Z.big[i], n)
            pkern <- dnorm(((zkk - new1)/band1), 0, 1)
            pkernel <- sum(pkern)/(n2 * band1)
            zkk <- rep(Z.big[i], n)
            pkern <- dnorm(((zkk - new2)/band2), 0, 1)
            pkernel2 <- sum(pkern)/(n2 * band2)
            if (i%%1000 == 0) 
                cat("number of estimated p-values:", i, fill = TRUE)
            c <- pkernel2/pkernel
            Ac.hat <- which(Tk.hat < c)
            Ac.hat.compl <- which(Tk.hat >= c)
            if (c <= min(Tk.hat)) 
                p.value[i] <- 0
            if (c >= max(Tk.hat)) 
                p.value[i] <- 1
            if (c > min(Tk.hat) && c < max(Tk.hat)) {
                integral <- (kern2[min(Ac.hat)] + sum(kern2[(min(Ac.hat) + 
                  1):(min(Ac.hat.compl) - 2)] * 2) + kern2[min(Ac.hat.compl) - 
                  1]) * integral.factor
                integral <- integral + (kern2[max(Ac.hat.compl) + 
                  1] + sum(kern2[(max(Ac.hat.compl) + 2):(max(Ac.hat) - 
                  1)] * 2) + kern2[max(Ac.hat)]) * integral.factor
                p.value[i] <- integral
            }
    if (is.na(p.value[i])) p.value[i] <- 0
        }
        cat(date(), fill = TRUE)
        return(p.value)
    }
}

