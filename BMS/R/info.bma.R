info.bma <-
function (object, ...) 
{
    bmao = object
    rm(object)
    foo = bmao$info
    iter = foo$iter
    burn = foo$burn
    timed = foo$timed
    models.visited = foo$models.visited
    corr.pmp = foo$corr.pmp
    K = foo$K
    N = foo$N
    msize = foo$msize
    cumsumweights = foo$cumsumweights
    if (is.element("mprior.info", names(bmao))) {
        prior = paste(bmao$mprior.info$mp.mode, "/", bmao$mprior.info$mp.msize)
    }
    else {
        if (is.element("theta", names(bmao$arguments)) && is.element("prior.msize", 
            names(bmao$arguments))) {
            if (!is.null(bmao$arguments$theta) & !is.null(bmao$arguments$prior.msize)) 
                prior = paste(bmao$arguments$theta, "/", bmao$arguments$prior.msize)
            else prior = NA
        }
        else {
            prior = paste(bmao$arguments$mprior, "/", bmao$arguments$mprior.size)
        }
    }
    gprior.info = bmao$gprior.info
    gprior.choice = gprior.info$gtype
    model.space = 2^K
    fraction.model = models.visited/model.space * 100
    fraction.topmodel = sum(bmao$topmod$ncount())/iter * 100
    if (gprior.info$gtype == "hyper") {
        gprior.choice = paste(gprior.choice, " (a=", 2 + signif(gprior.info$hyper.parameter - 
            2, digits = 4), ")", sep = "")
    }
    nr.reg = msize/cumsumweights
    info <- as.character(c(format(round(nr.reg, 4), nsmall = 4), 
        format(iter, nsmall = 0), format(burn, nsmall = 0), format(timed, 
            nsmall = 4), models.visited, format(model.space, 
            digits = 2), format(fraction.model, digits = 2), 
        format(fraction.topmodel, digits = 2), format(round(.cor.topmod(bmao$topmod), 
            4), nsmall = 4), format(N, nsmall = 4), prior, gprior.choice))
    names(info) <- c("Mean no. regressors", "Draws", "Burnins", 
        "Time", "No. models visited", "Modelspace 2^K", "% visited", 
        "% Topmodels", "Corr PMP", "No. Obs.", "Model Prior", 
        "g-Prior")
    if (gprior.info$return.g.stats) {
        gpriorav = gprior.info$shrinkage.moments[1]
        gstatsprint = paste("Av=", format(gpriorav, digits = 4), 
            sep = "")
        if (length(gprior.info$shrinkage.moments) > 1) {
            gpriorstdev = sqrt(gprior.info$shrinkage.moments[2] - 
                gprior.info$shrinkage.moments[1]^2)
            gstatsprint = paste(gstatsprint, ", Stdev=", format(gpriorstdev, 
                digits = 2), sep = "")
        }
        info <- c(info, gstatsprint)
        names(info)[13] <- "Shrinkage-Stats"
    }
    return(info)
}
