restricted <-
function (B, DD, yy, pp, lambda, smooth, nb, center, constmat, 
    types) 
{
    nterms = length(nb)
    m = length(yy)
    np = length(pp)
    lala <- matrix(lambda, nrow = nterms, ncol = 2, dimnames = list(1:nterms, 
        c("mean", "residual")))
    vector.a.ma.schall <- matrix(NA, nrow = sum(nb) + (1 * center), 
        ncol = np)
    if (smooth == "schall") {
        sch = schall(yy, B, 0.5, DD, nb, lala[, 1], constmat, 
            center, types)
        lala[, 1] = sch[[2]]
        mean.coefficients = sch[[1]]
        diag.hat = sch[[3]]
    }
    else if (smooth == "gcv") {
        acv.min = nlminb(start = lala[, 1], objective = acv, 
            yy = yy, B = B, quantile = 0.5, DD = DD, nb = nb, 
            constmat = constmat, lower = 0, upper = 10000)
        aa <- asyregpen.lsfit(yy, B, 0.5, abs(acv.min$par), DD, 
            nb, constmat)
        mean.coefficients <- aa$a
        lala[, 1] <- abs(acv.min$par)
        diag.hat = aa$diag.hat.ma
    }
    else if (smooth == "aic") {
        acv.min = nlminb(start = lala[, 1], objective = aicfun, 
            yy = yy, B = B, quantile = 0.5, DD = DD, nb = nb, 
            constmat = constmat, lower = 0, upper = 10000)
        aa <- asyregpen.lsfit(yy, B, 0.5, abs(acv.min$par), DD, 
            nb, constmat)
        mean.coefficients <- aa$a
        lala[, 1] <- abs(acv.min$par)
        diag.hat = aa$diag.hat.ma
    }
    else if (smooth == "bic") {
        acv.min = nlminb(start = lala[, 1], objective = bicfun, 
            yy = yy, B = B, quantile = 0.5, DD = DD, nb = nb, 
            constmat = constmat, lower = 0, upper = 10000)
        aa <- asyregpen.lsfit(yy, B, 0.5, abs(acv.min$par), DD, 
            nb, constmat)
        mean.coefficients <- aa$a
        lala[, 1] <- abs(acv.min$par)
        diag.hat = aa$diag.hat.ma
    }
    else if (smooth == "cvgrid") {
        lala[, 1] = cvgrid(yy, B, 0.5, DD, nb, constmat, types)
        aa <- asyregpen.lsfit(yy, B, 0.5, lala, DD, nb, constmat)
        mean.coefficients <- aa$a
        diag.hat = aa$diag.hat.ma
    }
    else if (smooth == "lcurve") {
        lala[, 1] = lcurve(yy, B, 0.5, DD, nb, constmat, types)
        aa <- asyregpen.lsfit(yy, B, 0.5, lala, DD, nb, constmat)
        mean.coefficients <- aa$a
        diag.hat = aa$diag.hat.ma
    }
    else {
        aa <- asyregpen.lsfit(yy, B, 0.5, lala[, 1], DD, nb, 
            constmat)
        mean.coefficients <- aa$a
        diag.hat = aa$diag.hat.ma
    }
    residuals = yy - B %*% mean.coefficients
    constmat[, ] = 0
    gg = asyregpen.lsfit(abs(residuals), B, 0.5, lala[, 1], DD, 
        nb, constmat)
    cc = NULL
    for (q in 1:np) {
        ca = asyregpen.lsfit(residuals, gg$fitted, pp[q], NULL, 
            matrix(0, nrow = 1, ncol = 1), 0, matrix(0, nrow = 1, 
                ncol = 1))
        cc[q] = ca$a
        vector.a.ma.schall[, q] = mean.coefficients + ca$a * 
            gg$a
    }
    diag.hat = matrix(diag.hat, nrow = length(diag.hat), ncol = np)
    return(list(vector.a.ma.schall, lala, diag.hat, mean.coefficients, 
        gg$a, cc))
}
