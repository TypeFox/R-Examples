logistol.int <- function (x, alpha = 0.05, P = 0.99, log.log = FALSE,
    side = 1) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (log.log) 
        x <- log(x)
    if (side == 2) {
        alpha <- alpha/2
	P <- (P + 1)/2
    }
    m.mom <- mean(x)
    s.mom <- sqrt(3 * (mean(x^2) - m.mom^2))/pi
    inits <- c(m.mom, s.mom)
    log.ll <- function(x, pars) sum(-dlogis(x, location = pars[1], 
        scale = pars[2], log = TRUE))
    out <- suppressWarnings(nlm(log.ll, p = inits, x = x, hessian = TRUE))
    out.est <- out$estimate
    m <- out.est[1]
    s <- out.est[2]
    inv.fish <- solve(out$hess)
    var.m <- inv.fish[1, 1]
    var.s <- inv.fish[2, 2]
    cov.ms <- inv.fish[1, 2]
    k.delta <- qlogis(P, scale = sqrt(3)/pi)
    t1 <- k.delta - cov.ms * qnorm(1 - alpha)^2
    t2 <- k.delta + cov.ms * qnorm(1 - alpha)^2
    u <- k.delta^2 - var.m * qnorm(1 - alpha)^2
    v <- 1 - var.s * qnorm(1 - alpha)^2
    k.lower <- (t1 + sqrt(t1^2 - u * v))/v
    k.upper <- (t2 + sqrt(t1^2 - u * v))/v
    L <- m - k.lower * s * pi/sqrt(3)
    U <- m + k.upper * s * pi/sqrt(3)
    if (log.log) {
        L <- exp(L)
        U <- exp(U)
    }
    if (side == 2) {
        alpha <- 2 * alpha
	P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, L, U))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "2-sided.lower", "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "1-sided.lower", "1-sided.upper")
    }
    temp
}

