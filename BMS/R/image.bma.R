image.bma <-
function (x, yprop2pip = FALSE, order.by.pip = TRUE, do.par = TRUE, 
    do.grid = TRUE, do.axis = TRUE, cex.axis = 1, ...) 
{
    dotargs = match.call(expand.dots = FALSE)$...
    ests = estimates.bma(x, exact = TRUE, order.by.pip = order.by.pip, 
        include.constant = FALSE)
    ests = ests[nrow(ests):1, ]
    pips = ests[, "PIP"]
    idx = ests[, "Idx"]
    pmp.res = pmp.bma(x, oldstyle = TRUE)
    pmps = pmp.res[, 1]
    normali_factor = sum(pmp.res[, 2])
    betasigns = beta.draws.bma(x)[idx, , drop = FALSE]
    betasigns = betasigns[as.logical(pips), ]
    betasigns = sign(betasigns)/2 + 0.5
    betasigns[betasigns == 0.5] = NA
    pips = pips[as.logical(pips)]
    if (yprop2pip) {
        pipbounds = (c(0, cumsum(pips)))
    }
    else {
        pipbounds = 0:length(pips)
        names(pipbounds) = c("", names(pips))
    }
    pmpbounds = (c(0, cumsum(pmps)))
    if (do.par) {
        oldmar = par()$mar
        spaceforyaxis = strwidth(names(pipbounds)[which.max(nchar(names(pipbounds)))], 
            units = "inches") * (par("mar")/par("mai"))[[2]]
        tempmar = oldmar
        tempmar[2] = min(spaceforyaxis + oldmar[2]/2, 0.5 * par("fin")[[1]] * 
            (par("mar")/par("mai"))[[2]])
        par(mar = tempmar)
    }
    dotargs = .adjustdots(dotargs, ylab = "", xlab = "Cumulative Model Probabilities", 
        col = c("tomato", "blue"), main = paste("Model Inclusion Based on Best ", 
            length(pmps), " Models"))
    dotargs$axes <- FALSE
    tbetasigns = t(betasigns)
    eval(as.call(c(list(as.name("image.default"), as.name("pmpbounds"), 
        as.name("pipbounds"), as.name("tbetasigns")), as.list(dotargs))))
    if (do.axis) {
        axis(1, at = pmpbounds, labels = round(normali_factor * 
            pmpbounds, 2), cex.axis = cex.axis)
        axis(2, at = pipbounds, labels = FALSE, line = FALSE)
        axis(2, at = pipbounds[-1] - diff(pipbounds)/2, labels = names(pipbounds[-1]), 
            tick = FALSE, las = 1, cex.axis = cex.axis)
    }
    if (do.grid) {
        abline(v = round(pmpbounds, 2), lty = "dotted", col = "grey")
        abline(h = round(pipbounds, 2), lty = "dotted", col = "grey")
    }
    if (do.par) {
        par(mar = oldmar)
    }
}
