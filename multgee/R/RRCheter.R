RRCheter <-
function (fmla, data, ncategories) 
{
    Consmat <- RCconstrains(ncategories, FALSE)
    dev <- stop.constrains <- Inf
    datax <- factor(data$x)
    datay <- factor(data$y)
    maxcategory <- nlevels(datax)
    noglm <- length(Consmat$nodf[Consmat$nodf < 2 * (maxcategory - 
        2)])
    fmla <- update(fmla, ~. - Mult(x, y) + Mult(z1, z2))
    for (i in 1:noglm) {
        data$z1 <- datax
        data$z2 <- datay
        levels(data$z1) <- pickcoefindz1 <- Consmat$parscores[i, 
            1:maxcategory]
        levels(data$z2) <- pickcoefindz2 <- Consmat$parscores[i, 
            -(1:maxcategory)]
        RRCmod <- suppressWarnings(gnm(fmla, data = data, family = poisson, 
            verbose = FALSE, model = FALSE))
        if (!is.null(RRCmod)) {
            if (deviance(RRCmod) < dev & RRCmod$conv) {
                scores <- as.numeric(coef(RRCmod)[pickCoef(RRCmod, 
                  "Mult")])
                pickcoefind <- unique(pickcoefindz1)
                scoresmu <- scores[pickcoefind][pickcoefindz1]
                mu <- normscores(scoresmu)
                if (all(diff(mu) >= 0) | all(diff(mu) <= 0)) {
                  scoresnu <- scores[-pickcoefind][pickcoefindz2]
                  nu <- normscores(scoresnu)
                  if (all(diff(nu) >= 0) | all(diff(nu) <= 0)) {
                    dev <- deviance(RRCmod)
                    stop.constrains <- Consmat$nodf[i]
                    LORterm <- c(tcrossprod(scoresmu, scoresnu))
                  }
                }
            }
        }
        if (stop.constrains < Consmat$nodf[i + 1]) 
            break
    }
    if (!is.finite(dev)) {
        fmla <- update(fmla, ~. - Mult(z1, z2) + x1:x2)
        for (i in (noglm + 1):length(Consmat$nodf)) {
            datax1 <- datax
            datay1 <- datay
            levels(datax1) <- pickcoefindz1 <- Consmat$parscores[i, 
                1:maxcategory]
            levels(datay1) <- pickcoefindz2 <- Consmat$parscores[i, 
                -(1:maxcategory)]
            data$x1 <- as.numeric(datax1)
            data$x2 <- as.numeric(datay1)
            RRCmod <- suppressWarnings(glm(fmla, data = data, 
                family = poisson))
            if (deviance(RRCmod) < dev & RRCmod$conv) {
                LORterm <- c(tcrossprod(pickcoefindz1, pickcoefindz2) * 
                  as.numeric(RRCmod$coef["x1:x2"]))
            }
        }
    }
    if (!is.finite(dev)) 
        LORterm <- rep(0, nlevels(datax)^2)
    LORterm
}

