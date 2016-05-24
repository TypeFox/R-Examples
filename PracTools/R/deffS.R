deffS <- function(p, w, y){
    if (any(w <= 0))
        warning("Some weights are less than or equal to 0.\n")
    if (any(p <= 0))
        warning("Some 1-draw selecti0on probabilities are less than or equal to 0.\n")

    n <- length(w)
    dK <- deffK(w)
    ybar <- weighted.mean(y,w)
    pbar <- weighted.mean(p,w)
    wbar <- weighted.mean(w,w)

    reg <- glm(y ~ p, weights = w)
    e <- reg$residuals
    A <- reg$coef[1]
    vy <- wtdvar(y,w)
    vp <- wtdvar(p,w)
    ve <- wtdvar(e,w)
    ve2 <- wtdvar(e^2,w)
    vw <- wtdvar(w,w)
    Nhat <- sum(w)
    ebar <- weighted.mean(e,w)
    e2bar <- weighted.mean(e^2,w)

    rho.yp <- sum(w*(y - ybar)*(p - pbar))/sum(w) / sqrt(vy * vp)
    rho.e2w <- sum(w*(e^2 - e2bar)*(w - wbar))/sum(w) / sqrt(ve2 * vw)
    rho.ew <- sum(w*(e - ebar)*(w - wbar))/sum(w) / sqrt(ve * vw)

    deffS <- A^2*(dK-1)/vy + dK*(1-rho.yp^2) + n*rho.e2w*sqrt(ve2*vw)/(Nhat*vy)
                + 2*n*A*rho.ew*sqrt(ve*vw)/(Nhat*vy)
    as.numeric(deffS)
}
