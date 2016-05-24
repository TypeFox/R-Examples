gdensity <-
function (x, n = 512, plot = TRUE, addons = "zles", addons.lwd = 1.5, 
    ...) 
{
    dsgivenykernel <- function(kpazvec, sf, N) {
        (kpazvec[[1]] - 2)/2 * (1 - sf)^((kpazvec[[1]] - 4)/2) * 
            (1 - sf * kpazvec[[2]])^(-(N - 1)/2)
    }
    if (!is.bma(x)) 
        stop("argument needs to an object of class 'bma'")
    if (!(x$gprior$gtype == "hyper")) 
        stop("g prior density makes only sense for hyper-g prior.")
    if (n < 2) 
        stop("n needs to be at least 2")
    n = floor(n)
    dotargs = match.call(expand.dots = FALSE)$...
    N = x$info$N
    K = x$info$K
    tm = x$topmod
    bools = tm$bool_binary()
    betas = tm$betas()
    betas2 = tm$betas2()
    smoments = tm$fixed_vector()
    yXdata = as.matrix(x$X.data)
    yXdata = yXdata - matrix(colMeans(yXdata), N, K + 1, byrow = TRUE)
    yty = c(crossprod(yXdata[, 1]))
    positions = lapply(lapply(as.list(as.data.frame(bools)), 
        as.logical), which)
    ymyvec = unlist(lapply(lapply(positions, .ols.terms2, yty = yty, 
        N = N, K = K, XtX.big = crossprod(yXdata[, -1]), Xty.big = c(crossprod(yXdata[, 
            -1], yXdata[, 1]))), function(x) x$full.results()$ymy))
    kvec = tm$kvec_raw()
    zvec = 1 - ymyvec/yty
    pmpexact = pmp.bma(x, oldstyle = TRUE)[, 1]
    f21a = x$gprior.info$hyper.parameter
    if (length(smoments) == 0) {
        lprob = x$gprior.info$lprobcalc
        smoments = sapply(lapply(as.list(as.data.frame(rbind(kvec, 
            ymyvec))), function(x) lprob$lprob.all(ymy = x[2], 
            k = x[1], bhat = numeric(x[1]), diag.inverse = rep(1, 
                x[1]))), "[[", "otherstats")
    }
    Es = c(crossprod(smoments[1, ], pmpexact))
    Es2 = c(crossprod(smoments[2, ], pmpexact))
    Esd = sqrt(Es2 - Es^2)
    nbsteps = n
    cutoff = max(0, Es - 5 * Esd)
    sdiff = (1 - cutoff)/(nbsteps + 1)
    s.seq = seq(sdiff + cutoff, cutoff + nbsteps * sdiff, sdiff)
    sdensl = lapply(as.list(as.data.frame(rbind(kvec + f21a, 
        zvec))), dsgivenykernel, sf = s.seq, N = N)
    intconsts = lapply(lapply(sdensl, sum), "*", sdiff)
    sdensvecs = mapply("/", sdensl, intconsts)
    sdens = sdensvecs %*% pmpexact
    reslist = list(x = s.seq, y = sdens, bw = NULL, n = n, call = sys.call(), 
        data.name = "Shrinkage", has.na = FALSE)
    class(reslist) = "density"
    if (!plot) {
        return(reslist)
    }
    dotargs = .adjustdots(dotargs, ylab = "Density", xlab = "Shrinkage factor", 
        main = "Posterior Density of the Shrinkage Factor", type = "l", 
        col = "steelblue4")
    eval(as.call(c(list(as.name("plot"), as.name("s.seq"), as.name("sdens")), 
        as.list(dotargs))))
    leg.col = numeric(0)
    leg.lty = numeric(0)
    leg.legend = character(0)
    if (any(grep("f", addons, ignore.case = TRUE))) {
        for (m in 1:length(pmpexact)) {
            Esm = smoments[1, m]
            if (as.logical(Esm)) {
                ixlower = max(sum(s.seq < Esm), 1)
                Esheight = (sdens[ixlower + 1] - sdens[ixlower]) * 
                  (Esm - s.seq[ixlower]) + sdens[ixlower]
                lines(x = rep(Esm, 2), y = c(0, Esheight), col = 8, 
                  lwd = addons.lwd)
            }
        }
        leg.col = c(leg.col, 8)
        leg.lty = c(leg.lty, 1)
        leg.legend = c(leg.legend, "EV Models")
    }
    if (any(grep("e", addons, ignore.case = FALSE))) {
        abline(v = Es, col = 2, lwd = addons.lwd)
        leg.col = c(leg.col, 2)
        leg.lty = c(leg.lty, 1)
        leg.legend = c(leg.legend, "EV")
    }
    if (any(grep("s", addons, ignore.case = FALSE))) {
        if (!(Es - 2 * Esd) < 0) 
            abline(v = Es - 2 * Esd, col = 2, lty = 2, lwd = addons.lwd)
        if (!(Es + 2 * Esd) > 1) 
            abline(v = Es + 2 * Esd, col = 2, lty = 2, lwd = addons.lwd)
        leg.col = c(leg.col, 2)
        leg.lty = c(leg.lty, 2)
        leg.legend = c(leg.legend, "2x SD")
    }
    if (any(grep("m", addons, ignore.case = TRUE))) {
        median_index = sum(cumsum(sdens) < sum(sdens)/2)
        abline(v = (s.seq[median_index] + s.seq[median_index + 
            1])/2, col = 3, lwd = addons.lwd)
        leg.col = c(leg.col, 3)
        leg.lty = c(leg.lty, 1)
        leg.legend = c(leg.legend, "Median")
    }
    if (any(grep("z", addons, ignore.case = TRUE))) {
        abline(h = 0, col = "gray", lwd = addons.lwd)
    }
    if (any(grep("E", addons, ignore.case = FALSE))) {
        if (all(x$gprior.info$shrinkage.moments == 0)) 
            warning("bma object needs to contain posterior g statistics - cf. argument 'g.stats' in 'help(bms)'")
        else {
            abline(v = x$gprior.info$shrinkage.moments[1], col = 4, 
                lwd = addons.lwd)
            leg.col = c(leg.col, 4)
            leg.lty = c(leg.lty, 1)
            leg.legend = c(leg.legend, "EV (MCMC)")
        }
    }
    if (any(grep("S", addons, ignore.case = FALSE))) {
        if (!all(x$gprior.info$shrinkage.moments == 0)) {
            ES = x$gprior.info$shrinkage.moments[1]
            SDs = sqrt(x$gprior.info$shrinkage.moments[2] - x$gprior.info$shrinkage.moments[1]^2)
            if (ES - 2 * SDs > 0) 
                abline(v = ES - 2 * SDs, col = 4, lty = 2, lwd = addons.lwd)
            if (ES + 2 * SDs < 1) 
                abline(v = ES + 2 * SDs, col = 4, lty = 2, lwd = addons.lwd)
            leg.col = c(leg.col, 4)
            leg.lty = c(leg.lty, 2)
            leg.legend = c(leg.legend, "2x SD (MCMC)")
        }
    }
    if (any(grep("l", addons, ignore.case = TRUE)) & (length(leg.col) > 
        0)) {
        leg.pos = "topleft"
        legend(x = leg.pos, lty = leg.lty, col = leg.col, legend = leg.legend, 
            box.lwd = 0, bty = "n")
    }
    return(invisible(reslist))
}
