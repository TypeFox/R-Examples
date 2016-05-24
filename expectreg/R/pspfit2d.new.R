pspfit2d.new <-
function (B, P, ps, yy, w, lambdax, lambdap, center) 
{
    x0 <- min(ps) - 0.001
    x1 <- max(ps) + 0.001
    B.size = length(ps)/nrow(B[[1]])
    B.deg = 1
    dx = (x1 - x0)/(B.size - 1)
    By = splineDesign(knots = seq(x0 - dx * B.deg, x1 + dx * 
        B.deg, by = dx), x = ps, ord = B.deg + 1)
    basisX = NULL
    penaltyX = NULL
    Bx = NULL
    for (k in 1:(length(yy)/nrow(B[[1]]))) {
        Bx = rbind(Bx, B[[1]])
    }
    if (center) {
        Bx = cbind(1, Bx)
        P[[1]] = rbind(0, cbind(0, P[[1]]))
    }
    nx <- ncol(Bx)
    ny <- ncol(By)
    B1 <- kronecker(Bx, t(rep(1, ny)))
    B2 <- kronecker(t(rep(1, nx)), By)
    basisX <- B1 * B2
    nb = ncol(basisX) - 1 * center
    Px <- kronecker(sqrt(lambdax[1]) * P[[1]], diag(ny))
    Dy <- diff(diag(ny), diff = 2)
    Py <- kronecker(diag(nx), sqrt(lambdap[1]) * Dy)
    penaltyX = rbind(Px, Py)
    if (length(B) > 1) 
        for (i in 2:length(B)) {
            Bx = NULL
            for (k in 1:(length(yy)/nrow(B[[i]]))) {
                Bx = rbind(Bx, B[[i]])
            }
            nx <- ncol(Bx)
            ny <- ncol(By)
            B1 <- kronecker(Bx, t(rep(1, ny)))
            B2 <- kronecker(t(rep(1, nx)), By)
            nb = c(nb, ncol(B1))
            basisX <- cbind(basisX, B1 * B2)
            Px <- kronecker(sqrt(lambdax[i]) * P[[i]], diag(ny))
            Dy <- diff(diag(ny), diff = 2)
            Py <- kronecker(diag(nx), sqrt(lambdap[i]) * Dy)
            Pxy = rbind(Px, Py)
            penaltyX = rbind(cbind(penaltyX, matrix(0, nrow = dim(penaltyX)[1], 
                ncol = dim(Pxy)[2])), cbind(matrix(0, nrow = dim(Pxy)[1], 
                ncol = dim(penaltyX)[2]), Pxy))
        }
    zx <- rep(0, nrow(penaltyX))
    zplus <- c(yy, zx)
    Bplus <- rbind(basisX, penaltyX)
    for (it in 1:10) {
        wplus <- c(w, 1 + zx)
        model <- lsfit(Bplus, zplus, intercept = FALSE, wt = wplus)
        pcoef = model$coeff
        pfit <- basisX %*% pcoef
        w0 <- w
        w <- ps * (yy > pfit) + (1 - ps) * (yy <= pfit)
        nw <- sum(w != w0)
        if (nw == 0) 
            break
        cat(it, nw, "\n")
    }
    pcoef = matrix(pcoef, nrow = B.size)
    if (center) {
        intercept = pcoef[, 1]
    }
    else intercept = 0
    my = nrow(B[[1]])
    np = length(ps)/my
    Z = list()
    Bg = basisX[, 1:nb[1]]
    ag = pcoef[1:nb[1]]
    Z[[1]] = Bg %*% ag
    dim(Z[[1]]) = c(my, np)
    if (length(B) > 1) 
        for (i in 2:length(B)) {
            Bg = NULL
            ag = NULL
            if (center) {
                Bx = NULL
                for (k in 1:(length(yy)/nrow(B[[i]]))) {
                  Bx = rbind(Bx, B[[i]])
                }
                Bx = cbind(1, Bx)
                nx <- ncol(Bx)
                ny <- ncol(By)
                B1 <- kronecker(Bx, t(rep(1, ny)))
                B2 <- kronecker(t(rep(1, nx)), By)
                Bg = B1 * B2
                ag = matrix(pcoef[(sum(nb[1:(i - 1)]) + 1):(sum(nb[1:i]))], 
                  nrow = B.size)
                ag = cbind(intercept, ag)
                ag = as.vector(ag)
            }
            else {
                Bg = basisX[, (sum(nb[1:(i - 1)]) + 1):(sum(nb[1:i]))]
                ag = pcoef[(sum(nb[1:(i - 1)]) + 1):(sum(nb[1:i]))]
            }
            Z[[i]] = Bg %*% ag
            dim(Z[[i]]) = c(my, np)
        }
    return(list(coef = t(pcoef), fit = pfit, hatma = hat(model$qr), 
        curves = Z, intercept = intercept))
}
