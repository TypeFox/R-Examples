deg <-
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
    kern1 <- kern(new1, z, lz)
    kern2 <- kern(new2, z, lz)
    if (max(kern2) <= max(kern1)) {
        print("No differentially genes.", quote = FALSE)
    }
    if (max(kern2) > max(kern1)) {
        cat("Determination of the rejection region", fill = TRUE)
        Tk.hat <- kern2[which(kern1 > 0)]/kern1[which(kern1 > 
            0)]
        single.alpha <- alpha / n
        integral.factor <- 0.005
        help <- 0
        c0 <- 0
        c1 <- max(Tk.hat)
        repeat {
            c <- (c1 + c0)/2
            Ac.hat <- which(Tk.hat < c)
            Ac.hat.compl <- which(Tk.hat >= c)
            integral <- 0
            integral <- (kern2[min(Ac.hat)] + sum(kern2[(min(Ac.hat) + 
                1):(min(Ac.hat.compl) - 2)] * 2) + kern2[min(Ac.hat.compl) - 
                1]) * integral.factor
            integral <- integral + (kern2[max(Ac.hat.compl) + 
                1] + sum(kern2[(max(Ac.hat.compl) + 2):(max(Ac.hat) - 
                1)] * 2) + kern2[max(Ac.hat)]) * integral.factor
            if (integral > single.alpha) 
                c1 <- c
            if (integral < single.alpha) 
                c0 <- c
            if (abs(integral - single.alpha) <= 1e-08) 
                break
            if (c <= 2e-06) 
                break
            if (help > 100) 
                break
            help <- help + 1
        }
        z <- z[which(kern1 > 0)]
        f <- which(Tk.hat >= c)
        region <- c(z[min(f)], z[max(f)])
        a <- which(Z.big < region[1])
        b <- which(Z.big > region[2])
        values <- sort(c(a, b))
        cat(date(), fill = TRUE)
        return(values)
    }
}

