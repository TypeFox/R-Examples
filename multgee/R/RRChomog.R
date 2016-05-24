RRChomog <-
function (fmla, data, ncategories) 
{
    Consmat <- RCconstrains(ncategories, TRUE)
    datax <- factor(data$x)
    datay <- factor(data$y)
    dev <- stop.constrains <- Inf
    noglm <- length(Consmat$nodf[Consmat$nodf < (nlevels(datax) - 
        2)])
    fmla <- update(fmla, ~. - MultHomog(x,y) + MultHomog(z1,z2))
    for (i in 1:noglm) {
        data$z1 <- datax
        data$z2 <- datay
        levels(data$z1) <- levels(data$z2) <- pickcoefind <- Consmat$parscores[i, 
            ]
        suppressWarnings(RRCmod <- gnm(fmla, family = poisson, 
            data = data, verbose = FALSE, model = FALSE))
        if (!is.null(RRCmod)) {
            if (deviance(RRCmod) < dev & RRCmod$conv) {
                pickcoef <- pickCoef(RRCmod, "MultHomog")
                scores <- as.numeric(coef(RRCmod)[pickcoef])[pickcoefind]
                mu <- normscores(scores)
                if (all(diff(mu) >= 0) | all(diff(mu) <= 0)) {
                  dev <- deviance(RRCmod)
                  LORterm <- c(tcrossprod(scores))
                  if (Consmat$nodf[i] < Consmat$nodf[i + 1]) 
                    break
                }
            }
        }
        if (stop.constrains < Consmat$nodf[i + 1]) 
            break
    }
    if (!is.finite(dev)) {
        fmla <- update(fmla, ~. - MultHomog(z1,z2) + 
            x1:x2)
        data$x1 <- as.numeric(datax1)
        data$x2 <- as.numeric(datay1)
        for (i in (noglm + 1):length(Consmat$nodf)) datax1 <- datax
        datay1 <- datay
        levels(datax1) <- levels(datay1) <- pickcoefind <- Consmat$parscores[i, 
            ]
        RRCmod <- suppressWarnings(glm(fmla, data = data, family = poisson, 
            verbose = FALSE, model = FALSE))
        if (deviance(RRCmod) < dev & RRCmod$conv) {
            LORterm <- c(tcrossprod(pickcoefind) * as.numeric(RRCmod$coef["x1:x2"]))
        }
    }
    if (!is.finite(dev)) 
        LORterm <- rep(0, nlevels(datax)^2)
    LORterm
}

