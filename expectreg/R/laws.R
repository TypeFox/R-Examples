laws <-
function (B, DD, yy, pp, lambda, smooth, nb, center, constmat, 
    types) 
{
    nterms = length(nb)
    m = length(yy)
    np = length(pp)
    if (length(lambda) < nterms) 
        lala = rep(lambda[1], nterms)
    else lala = lambda
    dummy.reg <- function(pp, lala, smooth, yy, B, DD, nb, nterms, 
        center) {
        cat("Expectile: ", pp, "\n")
        if (smooth == "schall") {
            sch = schall(yy, B, pp, DD, nb, lala, constmat, center, 
                types)
            lala = sch[[2]]
            vector.a.ma.schall = sch[[1]]
            diag.hat = sch[[3]]
        }
        else if (smooth == "gcv") {
            acv.min = nlminb(start = lala, objective = acv, yy = yy, 
                B = B, quantile = pp, DD = DD, nb = nb, constmat = constmat, 
                lower = 0, upper = 10000)
            aa <- asyregpen.lsfit(yy, B, pp, abs(acv.min$par), 
                DD, nb, constmat)
            vector.a.ma.schall <- aa$a
            lala <- abs(acv.min$par)
            diag.hat = aa$diag.hat.ma
        }
        else if (smooth == "aic") {
            acv.min = nlminb(start = lala, objective = aicfun, 
                yy = yy, B = B, quantile = pp, DD = DD, nb = nb, 
                constmat = constmat, lower = 0, upper = 10000)
            aa <- asyregpen.lsfit(yy, B, pp, abs(acv.min$par), 
                DD, nb, constmat)
            vector.a.ma.schall <- aa$a
            lala <- abs(acv.min$par)
            diag.hat = aa$diag.hat.ma
        }
        else if (smooth == "bic") {
            acv.min = nlminb(start = lala, objective = bicfun, 
                yy = yy, B = B, quantile = pp, DD = DD, nb = nb, 
                constmat = constmat, lower = 0, upper = 10000)
            aa <- asyregpen.lsfit(yy, B, pp, abs(acv.min$par), 
                DD, nb, constmat)
            vector.a.ma.schall <- aa$a
            lala <- abs(acv.min$par)
            diag.hat = aa$diag.hat.ma
        }
        else if (smooth == "cvgrid") {
            lala = cvgrid(yy, B, pp, DD, nb, constmat, types)
            aa <- asyregpen.lsfit(yy, B, pp, lala, DD, nb, constmat)
            vector.a.ma.schall <- aa$a
            diag.hat = aa$diag.hat.ma
        }
        else if (smooth == "lcurve") {
            lala = lcurve(yy, B, pp, DD, nb, constmat, types)
            aa <- asyregpen.lsfit(yy, B, pp, lala, DD, nb, constmat)
            vector.a.ma.schall <- aa$a
            diag.hat = aa$diag.hat.ma
        }
        else {
            aa <- asyregpen.lsfit(yy, B, pp, lala, DD, nb, constmat)
            vector.a.ma.schall <- aa$a
            diag.hat = aa$diag.hat.ma
        }
        list(vector.a.ma.schall, lala, diag.hat)
    }
    if (.Platform$OS.type == "unix") 
        coef.vector = mclapply(pp, function(pp) dummy.reg(pp, 
            lala, smooth, yy, B, DD, nb, nterms, center), mc.cores = max(1, 
            min(detectCores() - 1, 2)))
    else if (.Platform$OS.type == "windows") 
        coef.vector = mclapply(pp, function(pp) dummy.reg(pp, 
            lala, smooth, yy, B, DD, nb, nterms, center), mc.cores = 1)
    lala <- matrix(lambda, nrow = nterms, ncol = np)
    vector.a.ma.schall <- matrix(NA, nrow = sum(nb) + (1 * center), 
        ncol = np)
    diag.hat = matrix(NA, nrow = m, ncol = np)
    for (i in 1:np) {
        vector.a.ma.schall[, i] = coef.vector[[i]][[1]]
        lala[, i] = coef.vector[[i]][[2]]
        diag.hat[, i] = coef.vector[[i]][[3]]
    }
    return(list(vector.a.ma.schall, lala, diag.hat))
}
