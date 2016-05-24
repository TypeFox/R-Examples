density.bma <-
function (x, reg = NULL, addons = "lemsz", std.coefs = FALSE, 
    n = 300, plot = TRUE, hnbsteps = 30, addons.lwd = 1.5, ...) 
{
    dtcm = function(x, df, ncp, varp) {
        sqvarp = sqrt(varp)
        dt((x - ncp)/sqvarp, df = df)/sqvarp
    }
    dsgivenykernel <- function(sf, kpa, N, z) {
        (kpa - 2)/2 * (1 - sf)^((kpa - 4)/2) * (1 - sf * z)^(-(N - 
            1)/2)
    }
    dotargs = match.call(expand.dots = FALSE)$...
    bmao = x
    if (!is.bma(bmao)) 
        stop("Argument bmao needs to be a bma object")
    if (hnbsteps%%2) 
        stop("Argument nbsteps needs to be an even integer")
    nbsteps = max(hnbsteps, 2)
    n = max(ceiling(n), 1)
    N = bmao$info$N
    K = bmao$info$K
    if (is.null(reg)) 
        reg = 1:K
    nameix = 1:K
    names(nameix) = bmao$reg.names
    reg = nameix[reg]
    ishyper = (bmao$gprior$gtype == "hyper")
    tm = bmao$topmod
    bools = (tm$bool_binary())
    betas = tm$betas()
    betas2 = tm$betas2()
    if (std.coefs) {
        sddata = apply(as.matrix(bmao$X.data), 2, stats::sd)
        betas = diag(sddata[-1]) %*% betas/sddata[1]
        betas2 = diag(sddata[-1]^2) %*% betas2/sddata[1]^2
    }
    sigmadiag = (betas2 - betas^2) * (N - 3)/(N - 1)
    pmps = pmp.bma(bmao$topmod, oldstyle = TRUE)[, 1]
    pips = c(tcrossprod(bools, t(pmps)))
    Eb1 = c(tcrossprod(betas, t(pmps)))/pips
    Ebsd = sqrt(c(tcrossprod(betas2, t(pmps)))/pips - Eb1^2)
    Ebsd[is.nan(Ebsd)] = 0
    Eb1[is.nan(Eb1)] = 0
    Eball = cbind(Eb1, Ebsd)
    if ((any(grep("E", addons, ignore.case = FALSE))) | (any(grep("S", 
        addons, ignore.case = FALSE)))) {
        Eb1.mcmc = bmao$info$b1mo/bmao$info$inccount
        Ebsd.mcmc = sqrt(bmao$info$b2mo/bmao$info$inccount - 
            Eb1.mcmc^2)
        if (std.coefs) {
            sddata = apply(as.matrix(bmao$X.data), 2, stats::sd)
            Eb1.mcmc = Eb1.mcmc * sddata[-1]/sddata[1]
            Ebsd.mcmc = Ebsd.mcmc * sddata[-1]/sddata[1]
        }
    }
    if (ishyper) {
        yXdata = as.matrix(bmao$X.data)
        yXdata = yXdata - matrix(colMeans(yXdata), N, K + 1, 
            byrow = TRUE)
        if (std.coefs) 
            yXdata = yXdata %*% diag(1/sddata)
        yty = c(crossprod(yXdata[, 1]))
        positions = lapply(lapply(as.list(as.data.frame(bools)), 
            as.logical), which)
        olsmodels = lapply(lapply(positions, .ols.terms2, yty = yty, 
            N = N, K = K, XtX.big = crossprod(yXdata[, -1]), 
            Xty.big = c(crossprod(yXdata[, -1], yXdata[, 1]))), 
            function(x) x$full.results())
        f21a = bmao$gprior.info$hyper.parameter
    }
    plotndens <- function(ix, doplot = FALSE) {
        sss = function(lbound, uboundp1, nbsteps) {
            s.seq = seq(lbound, uboundp1, (uboundp1 - lbound)/nbsteps)[-nbsteps]
            tmat = sapply(as.list(s.seq), function(ss) {
                dtcm(seqs, N - 1, ss * bhati, invdiagi * ss * 
                  (1 - ss * z)/(N - 1) * yty)
            })
            smat = sapply(as.list(s.seq), dsgivenykernel, kpa = k + 
                f21a, N = N, z = z)
            if (any(is.infinite(smat))) 
                smat[is.infinite(smat)] = 0
            intconst = (4 * sum(smat[c(FALSE, TRUE)]) + 2 * sum(smat[c(TRUE, 
                FALSE)]) - 3 * smat[nbsteps] - smat[1]) * (s.seq[nbsteps] - 
                s.seq[1])/nbsteps/3
            return(list(dv = c(4 * tmat[, c(FALSE, TRUE)] %*% 
                smat[c(FALSE, TRUE)] + 2 * tmat[, c(TRUE, FALSE)] %*% 
                smat[c(TRUE, FALSE)] - 3 * tmat[, nbsteps] * 
                smat[nbsteps] - tmat[, 1] * smat[1]) * (s.seq[nbsteps] - 
                s.seq[1])/nbsteps/3, ic = intconst))
        }
        if (pips[ix] == 0) {
            reslist = list(x = numeric(n), y = numeric(n), n = n, 
                call = sys.call(), data.name = names(nameix)[ix], 
                has.na = FALSE)
            class(reslist) = c("density", "coef.density")
            return(reslist)
        }
        lbound = min(betas[ix, as.logical(bools[ix, ])]) - 3 * 
            Eball[ix, 2]
        ubound = max(betas[ix, as.logical(bools[ix, ])]) + 3 * 
            Eball[ix, 2]
        seqs = seq(lbound, ubound, (ubound - lbound)/(n - 1))
        densvec = numeric(length(seqs))
        for (m in 1:length(pmps)) {
            if (bools[ix, m]) {
                if (ishyper) {
                  ixadj = sum(bools[1:ix, m])
                  bhati = olsmodels[[m]]$bhat[[ixadj]]
                  invdiagi = olsmodels[[m]]$diag.inverse[[ixadj]]
                  k = sum(bools[, m])
                  Esf = betas[ix, m]/bhati
                  z = 1 - olsmodels[[m]]$ymy/yty
                  midpoint = 1 - (1 - Esf) * 4
                  if (midpoint < 0.5) {
                    dvl = sss(1e-04, 0.9999999, nbsteps * 2)
                    addvec = dvl$dv/dvl$ic
                  }
                  else {
                    dvl1 = sss(1e-04, midpoint, nbsteps)
                    dvl2 = sss(midpoint, 1, nbsteps)
                    addvec = (dvl1$dv + dvl2$dv)/(dvl1$ic + dvl2$ic)
                  }
                }
                else {
                  addvec = dtcm(seqs, N - 1, betas[ix, m], sigmadiag[ix, 
                    m])
                }
                densvec = densvec + pmps[m] * addvec
            }
        }
        reslist = list(x = seqs, y = densvec, bw = NULL, n = n, 
            call = sys.call(), data.name = names(nameix)[ix], 
            has.na = FALSE)
        class(reslist) = "density"
        if (!doplot) {
            return(reslist)
        }
        main_default = paste("Marginal Density:", names(nameix)[ix], 
            "(PIP", round(c(crossprod(pmps, bools[ix, ])) * 100, 
                2), "%)")
        if (any(grep("p", addons, ignore.case = TRUE))) {
            decr = 0.12
            parplt = par()$plt
            parplt_temp = parplt
            parplt_temp[4] = (1 - decr) * parplt[4] + decr * 
                parplt[3]
            par(plt = parplt_temp)
            main_temp = main_default
            main_default = NULL
        }
        dotargs = .adjustdots(dotargs, type = "l", col = "steelblue4", 
            main = main_default, xlab = if (std.coefs) 
                "Standardized Coefficient"
            else "Coefficient", ylab = "Density")
        eval(as.call(c(list(as.name("plot"), x = as.name("seqs"), 
            y = as.name("densvec")), as.list(dotargs))))
        leg.col = numeric(0)
        leg.lty = numeric(0)
        leg.legend = character(0)
        if (any(grep("g", addons, ignore.case = TRUE))) {
            grid()
        }
        if (any(grep("b", addons, ignore.case = TRUE))) {
            for (m in 1:length(pmps)) {
                Ebm = betas[ix, m]
                if (as.logical(Ebm)) {
                  Ebheight = min(densvec[max(sum(seqs < Ebm), 
                    1)], densvec[sum(seqs < Ebm) + 1])
                  lines(x = rep(Ebm, 2), y = c(0, Ebheight), 
                    col = 8)
                }
            }
            leg.col = c(leg.col, 8)
            leg.lty = c(leg.lty, 1)
            leg.legend = c(leg.legend, "EV Models")
        }
        if (any(grep("e", addons, ignore.case = FALSE))) {
            abline(v = Eball[ix, 1], col = 2, lwd = addons.lwd)
            leg.col = c(leg.col, 2)
            leg.lty = c(leg.lty, 1)
            leg.legend = c(leg.legend, "Cond. EV")
        }
        if (any(grep("s", addons, ignore.case = FALSE))) {
            abline(v = Eball[ix, 1] - 2 * Eball[ix, 2], col = 2, 
                lty = 2, lwd = addons.lwd)
            abline(v = Eball[ix, 1] + 2 * Eball[ix, 2], col = 2, 
                lty = 2, lwd = addons.lwd)
            leg.col = c(leg.col, 2)
            leg.lty = c(leg.lty, 2)
            leg.legend = c(leg.legend, "2x Cond. SD")
        }
        if (any(grep("m", addons, ignore.case = TRUE))) {
            median_index = sum(cumsum(densvec) < sum(densvec)/2)
            abline(v = (seqs[median_index] + seqs[median_index + 
                1])/2, col = 3, lwd = addons.lwd)
            leg.col = c(leg.col, 3)
            leg.lty = c(leg.lty, 1)
            leg.legend = c(leg.legend, "Median")
        }
        if (any(grep("z", addons, ignore.case = TRUE))) {
            abline(h = 0, col = "gray", lwd = addons.lwd)
        }
        if (any(grep("E", addons, ignore.case = FALSE))) {
            abline(v = Eb1.mcmc[ix], col = 4, lwd = addons.lwd)
            leg.col = c(leg.col, 4)
            leg.lty = c(leg.lty, 1)
            leg.legend = c(leg.legend, "Cond. EV (MCMC)")
        }
        if (any(grep("S", addons, ignore.case = FALSE))) {
            abline(v = Eb1.mcmc[ix] - 2 * Ebsd.mcmc[ix], col = 4, 
                lty = 2, lwd = addons.lwd)
            abline(v = Eb1.mcmc[ix] + 2 * Ebsd.mcmc[ix], col = 4, 
                lty = 2, lwd = addons.lwd)
            leg.col = c(leg.col, 4)
            leg.lty = c(leg.lty, 2)
            leg.legend = c(leg.legend, "2x SD (MCMC)")
        }
        if (any(grep("l", addons, ignore.case = TRUE)) & (length(leg.col) > 
            0)) {
            leg.pos = "topright"
            if (Eball[ix, 1] > seqs[floor(n/2)]) 
                leg.pos = "topleft"
            legend(x = leg.pos, lty = leg.lty, col = leg.col, 
                legend = leg.legend, box.lwd = 0, bty = "n", 
                lwd = addons.lwd)
        }
        if (any(grep("p", addons, ignore.case = TRUE))) {
            pusr = par()$usr
            rect(pusr[1], pusr[4] * (1 + decr * 0.2), pusr[2], 
                pusr[4] * (1 + decr), xpd = TRUE, col = 8)
            rect(pusr[1], pusr[4] * (1 + decr * 0.2), pips[ix] * 
                pusr[2] + (1 - pips[ix]) * pusr[1], pusr[4] * 
                (1 + decr), xpd = TRUE, col = 9)
            mtext("PIP:", side = 2, las = 2, line = 1, at = pusr[4] * 
                (1 + decr * 0.6))
            par(plt = parplt)
            title(main_temp)
        }
        return(reslist)
    }
    densres = list()
    oldask = par()$ask
    plots = 0
    for (vbl in 1:length(reg)) {
        doplot = (if (as.logical(pips[reg[vbl]])) 
            plot
        else FALSE)
        plots = plots + doplot
        if (plots == 2) {
            par(ask = TRUE)
        }
        densres[[nameix[vbl]]] = plotndens(reg[vbl], doplot)
        densres[[nameix[vbl]]]$call = sys.call()
    }
    par(ask = oldask)
    if (length(densres) == 1) 
        densres = densres[[1]]
    else class(densres) = c("coef.density", class(densres))
    if (!plot) 
        return(densres)
    if (plot & (plots == 0)) {
        warning("No plot produced as PIPs of provided variables are zero under 'exact' estimation.")
    }
    return(invisible(densres))
}
