deffH <- function(w, y, x){
    if (any(w <= 0))
        warning("Some weights are less than or equal to 0.\n")

    n <- length(w)
    dK <- deffK(w)
    wbar <- weighted.mean(w,w)

    reg <- glm(y ~ x, weights = w)
    e <- reg$residuals
    u <- reg$coef[1] + e
    A <- reg$coef[1]

    vy <- wtdvar(y,w)
    vu <- wtdvar(u,w)
    vu2 <- wtdvar(u^2,w)
    vw <- wtdvar(w,w)
    Nhat <- sum(w)
    ubar <- weighted.mean(u,w)
    u2bar <- weighted.mean(u^2,w)

    rho.u2w <- sum(w*(u^2 - u2bar)*(w - wbar))/sum(w) / sqrt(vu2 * vw)
    rho.uw <- sum(w*(u - ubar)*(w - wbar))/sum(w) / sqrt(vu * vw)

    deffH <- dK * vu/vy + n*sqrt(vw)* (rho.u2w*sqrt(vu) - 2*A*rho.uw*sqrt(vu)) / (Nhat*vy)
    as.numeric(deffH)
}
