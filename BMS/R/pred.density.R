pred.density <-
function (object, newdata = NULL, n = 300, hnbsteps = 30, ...) 
{
    dtcm = function(x, df, ncp, varp) {
        sqvarp = sqrt(varp)
        dt((x - ncp)/sqvarp, df = df)/sqvarp
    }
    dsgivenykernel <- function(sf, kpa, N, z) {
        (kpa - 2)/2 * (1 - sf)^((kpa - 4)/2) * (1 - sf * z)^(-(N - 
            1)/2)
    }
    nbsteps = max(hnbsteps, 2)
    n = max(ceiling(n), 1)
    is.hyper = (object$gprior.info$gtype == "hyper")
    if (is.hyper) 
        f21a = object$gprior.info$hyper.parameter
    if (is.bma(object)) {
        K = object$info$K
        N = object$info$N
        yXdata = as.matrix(object$X.data)
        tmo <- object$topmod
    }
    else if (is(object, "zlm")) {
        yXdata = as.matrix(object$model)
        K = ncol(yXdata) - 1
        N = nrow(yXdata)
        tmo <- topmod(1, nmaxregressors = K, bbeta = TRUE, liks = object$marg.lik, 
            ncounts = 1, modelbinaries = matrix(rep(1, K), K, 
                1), betas = matrix(as.vector(object$coefficients[-1]), 
                K), betas2 = matrix(as.vector(object$coef2moments[-1]), 
                K))
    }
    else stop("argument 'object' requires class 'bma' or 'zlm'")
    rm(object)
    if (missing(newdata)) {
        stop("You must provide the argument newdata")
    }
    else {
        newX = as.matrix(newdata)
        if (!is.numeric(newX)) 
            stop("newdata must be numeric!")
        if (is.vector(newdata)) 
            newX = matrix(newdata, 1)
        if (ncol(newX) != K) {
            if (ncol(newX) == K + 1) {
                newX = newX[, -1, drop = FALSE]
            }
            else {
                stop("newdata must be a matrix or data.frame with ", 
                  K, " columns.")
            }
        }
        orinames = colnames(yXdata[, -1, drop = FALSE])
        if (!is.null(colnames(newX)) && !is.null(orinames)) {
            if (all(orinames %in% colnames(newX)) && !all(orinames == 
                colnames(newX))) {
                warning("argument newdata had to be reordered according to its column names. Consider submitting the columns of newdata in the right order.")
                newX = newX[, orinames, drop = FALSE]
            }
        }
    }
    if (!is.null(rownames(newX))) {
        newXnames = rownames(newX)
    }
    else {
        newXnames = as.character(1:nrow(newX))
    }
    rnew = nrow(newX)
    y.mean = mean(yXdata[, 1])
    y <- yXdata[, 1] - matrix(y.mean, N, 1, byrow = TRUE)
    X <- yXdata[, -1, drop = FALSE] - matrix(colMeans(yXdata[, 
        -1, drop = FALSE]), N, K, byrow = TRUE)
    XtX.big = crossprod(X)
    Xty.big = as.vector(crossprod(X, y))
    yty = crossprod(y)[[1]]
    newXdm = newX - matrix(colMeans(yXdata[, -1, drop = FALSE]), 
        rnew, K, byrow = TRUE)
    hexobject <- .hexcode.binvec.convert(K)
    make_xfxxxf = function(hex) {
        syminv <- function(symmat, ndim = ncol(symmat)) {
            if (!is.matrix(symmat)) {
                symmat = as.matrix(symmat)
            }
            if (dim(symmat)[[1]] == 0) 
                return(matrix(numeric(0), 0, 0))
            return(chol2inv(chol(symmat), size = ndim))
        }
        boolvec = as.logical(hexobject$as.binvec(hex))
        if (!any(boolvec)) 
            return(c(numeric(rnew), numeric(rnew), Inf, Inf, 
                0))
        newXsub = newXdm[, boolvec, drop = FALSE]
        xtxinv = syminv(XtX.big[boolvec, boolvec, drop = FALSE])
        xty = Xty.big[boolvec]
        betas = as.vector(crossprod(xtxinv, xty), mode = "numeric")
        r2 = crossprod(xty, betas)[[1]]/yty
        xtxinv_xf = tcrossprod(xtxinv, newXsub)
        xf_xx_xf = unlist(lapply(1:nrow(newXsub), function(x) {
            crossprod(newXsub[x, ], xtxinv_xf[, x])[[1L]]
        }))
        xf_bhat = as.vector(newXsub %*% betas)
        return(c(xf_xx_xf, xf_bhat, xtxinv[[1L]], betas[[1L]], 
            r2))
    }
    pmps = pmp.bma(tmo, oldstyle = TRUE)[, 1, drop = TRUE]
    bools = tmo$bool()
    nmodel = length(bools)
    linvres = lapply(bools, make_xfxxxf)
    mat_xfxxxf = array(unlist(lapply(linvres, "[", 1:rnew)), 
        dim = c(rnew, nmodel))
    mat_xfbhat = array(unlist(lapply(linvres, "[", rnew + (1:rnew))), 
        dim = c(rnew, nmodel))
    xtxinv_elem1 = unlist(lapply(linvres, "[[", rnew * 2 + 1))
    betahat_elem1 = unlist(lapply(linvres, "[[", rnew * 2 + 2))
    r2 = unlist(lapply(linvres, "[[", rnew * 2 + 3))
    kvec = tmo$kvec_raw()
    kvec_cs = c(1, cumsum(kvec) + 1)
    kvec_cs = kvec_cs[-length(kvec_cs)]
    firstbetas = tmo$betas_raw()[kvec_cs]
    firstbetas2 = tmo$betas2_raw()[kvec_cs]
    Es = firstbetas/betahat_elem1
    varmult = (firstbetas2 - firstbetas^2)/xtxinv_elem1
    if (is.hyper) {
        first_factor = yty/(N - 3) * (N + 1)/N * (1 + 2/(N - 
            kvec - f21a - 1) - r2 * Es)
    }
    else {
        first_factor = yty/(N - 3) * (1 - Es * r2) * (N + 1)/N
    }
    Sigmas = (matrix((N - 3)/(N - 1) * first_factor, rnew, nmodel, 
        byrow = TRUE) + t(t(mat_xfxxxf) * ((N - 3)/(N - 1) * 
        varmult)))
    Evals_minusy = t(t(mat_xfbhat) * Es)
    Eyf = as.vector(Evals_minusy %*% pmps + y.mean)
    Varyf = as.vector(Sigmas %*% pmps) * (N - 1)/(N - 3)
    premultfactor = yty/(N - 1)
    interceptfactor = (N + 1)/N
    calcdensvec = function(xf_index, seqy, m_index) {
        sss = function(lbound, uboundp1, nbsteps, seqs, xf.index) {
            s.seq = seq(lbound, uboundp1, (uboundp1 - lbound)/nbsteps)[-nbsteps]
            tmat = array(unlist(lapply(as.list(s.seq), function(ss) {
                dtcm(seqs, N - 1, y.mean + ss * myev, premultfactor * 
                  (1 - ss * myr2) * (interceptfactor + ss * myxfxxxf))
            })), dim = c(length(seqs), nbsteps))
            smat = sapply(as.list(s.seq), dsgivenykernel, kpa = myk + 
                f21a, N = N, z = myr2)
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
        if (any(is.na(newX[xf_index, ]))) {
            densvec = numeric(0)
        }
        if (is.hyper) {
            myev = mat_xfbhat[xf_index, m_index]
            myxfxxxf = mat_xfxxxf[xf_index, m_index]
            myk = kvec[[m_index]]
            myr2 = r2[[m_index]]
            midpoint = 1 - (1 - Es[[m_index]]) * 4
            if (midpoint < 0.5) {
                dvl = sss(1e-04, 0.9999999, nbsteps * 2, seqy, 
                  xf_index)
                densvec = dvl$dv/dvl$ic
            }
            else {
                dvl1 = sss(1e-04, midpoint, nbsteps, seqy, xf_index)
                dvl2 = sss(midpoint, 1, nbsteps, seqy, xf_index)
                densvec = (dvl1$dv + dvl2$dv)/(dvl1$ic + dvl2$ic)
            }
        }
        else {
            densvec = dtcm(seqy, N - 1, Evals_minusy[xf_index, 
                m_index] + y.mean, Sigmas[xf_index, m_index])
        }
        return(densvec)
    }
    dens_yf = function(yfr, xf_indices = NULL) {
        if (is.null(xf_indices)) 
            xf_indices = seq_len(rnew)
        yfdens = array(NA, dim = dim(yfr))
        for (myxf in 1:length(xf_indices)) {
            allm_dens = sapply(seq_len(nmodel), function(x) calcdensvec(xf_indices[[myxf]], 
                yfr[myxf, ], x))
            yfdens[myxf, ] = as.vector(allm_dens %*% pmps)
        }
        yfdens[!is.finite(yfdens)] = NA
        if (ncol(yfdens) == 1) 
            dim(yfdens) <- NULL
        return(yfdens)
    }
    emptydens = list(x = numeric(0), y = numeric(0), bw = NULL, 
        n = 0, has.na = TRUE)
    class(emptydens) = "density"
    dlist = lapply(vector("list", nrow(newX)), function(x) emptydens)
    densities_calculated <- FALSE
    calc_alldens = function() {
        if (densities_calculated) 
            return(NULL)
        for (xf.index in 1:rnew) {
            if (!any(is.na(newX[xf.index, ]))) {
                lbound = Eyf[[xf.index]] - sqrt(Varyf[[xf.index]]) * 
                  4
                ubound = Eyf[[xf.index]] + sqrt(Varyf[[xf.index]]) * 
                  4
                seqs = seq(lbound, ubound, (ubound - lbound)/(n - 
                  1))
                allm_dens = sapply(seq_len(nmodel), function(x) calcdensvec(xf.index, 
                  seqs, x))
                myy = as.vector(tcrossprod(t(as.matrix(pmps)), 
                  allm_dens))
                mydens = list(x = seqs, y = myy, bw = NULL, n = n, 
                  has.na = FALSE)
                class(mydens) = "density"
                dlist[[xf.index]] <<- mydens
            }
        }
        densities_calculated <<- TRUE
    }
    consistent.yf = function(yf, xf.indices = NULL) {
        xf_series = seq_len(rnew)
        wasnull = FALSE
        if (is.null(xf.indices)) {
            wasnull = TRUE
            xf.indices = xf_series
        }
        else {
            if (!all(xf.indices %in% xf_series)) 
                stop(paste("predict_index needs to be an integer between 1 and ", 
                  rnew, "!", sep = ""))
        }
        if (!is.numeric(yf)) 
            stop("realized.y must be a numeric matrix or vector!")
        if (!is.matrix(yf)) 
            yf <- as.matrix(yf)
        if ((length(xf.indices) == 1) & (nrow(yf) > 1) & (ncol(yf) == 
            1)) 
            yf <- t(yf)
        if (nrow(newX[xf.indices, , drop = FALSE]) != nrow(yf)) {
            if (wasnull) 
                stop(paste("realized.y must have", rnew, "elements/rows corresponding to newdata"))
            else stop("The number of rows/elements in realized.y must have the same length as predict_index!")
        }
        return(yf)
    }
    consistent_predict_index = function(pix) {
        if (is.character(pix)) {
            if (all(pix %in% newXnames)) {
                return(match(pix, newXnames))
            }
            else {
                stop("Forecast IDs provided in predict_index do not conform to rownames of predicted data")
            }
        }
        else return(pix)
    }
    plot.preddens = function(xf.index = 1, addons = "eslz", yf.addons = NULL, 
        predict_index = NULL, addons.lwd = 1.5, ...) {
        dotargs = match.call(expand.dots = FALSE)$...
        if (rnew > 1) {
            main_default <- paste("Predictive Density Obs ", 
                newXnames[[xf.index]], " (", nmodel, " Models)", 
                sep = "")
        }
        else {
            main_default <- paste("Predictive Density", " (", 
                nmodel, " Models)", sep = "")
        }
        dotargs = .adjustdots(dotargs, xlab = "Response variable", 
            main = main_default, col = 4, zero.line = FALSE)
        thingy = dlist[[xf.index]]
        eval(as.call(c(list(as.name("plot"), as.name("thingy")), 
            as.list(dotargs))))
        leg.col = numeric(0)
        leg.lty = numeric(0)
        leg.legend = character(0)
        if (any(grep("g", addons, ignore.case = TRUE))) {
            grid()
        }
        if (any(grep("e", addons, ignore.case = FALSE))) {
            abline(v = fit[[xf.index]], col = 2, lwd = addons.lwd)
            leg.col = c(leg.col, 2)
            leg.lty = c(leg.lty, 1)
            leg.legend = c(leg.legend, "Exp. Value")
        }
        if (any(grep("s", addons, ignore.case = FALSE))) {
            abline(v = fit[[xf.index]] - 2 * stderrs[[xf.index]], 
                col = 2, lty = 2, lwd = addons.lwd)
            abline(v = fit[[xf.index]] + 2 * stderrs[[xf.index]], 
                col = 2, lty = 2, lwd = addons.lwd)
            leg.col = c(leg.col, 2)
            leg.lty = c(leg.lty, 2)
            leg.legend = c(leg.legend, "2x Std.Errs")
        }
        if (any(grep("z", addons, ignore.case = TRUE))) {
            abline(h = 0, col = "gray", lwd = addons.lwd)
        }
        if (!is.null(yf.addons) && is.numeric(yf.addons)) {
            yfs = as.vector(yf.addons)
            if (!is.na(yfs[[xf.index]])) {
                abline(v = yfs[[xf.index]], col = 1, lwd = addons.lwd, 
                  lty = 2)
                leg.col = c(leg.col, 1)
                leg.lty = c(leg.lty, 2)
                leg.legend = c(leg.legend, "Realized y")
            }
            else warning("yf.addons must be a vector with the same number of elements as rows in newdata!")
        }
        if (any(grep("l", addons, ignore.case = TRUE)) & (length(leg.col) > 
            0)) {
            legend(x = "topright", lty = leg.lty, col = leg.col, 
                legend = leg.legend, box.lwd = 0, bty = "n", 
                lwd = addons.lwd)
        }
    }
    fit = Eyf
    names(fit) = newXnames
    stderrs = sqrt(Varyf)
    names(stderrs) = newXnames
    reslist = list()
    reslist$densities = function() {
        calc_alldens()
        return(dlist)
    }
    reslist$fit = fit
    reslist$std.err = stderrs
    reslist$dyf = function(realized.y, predict_index = NULL) {
        predict_index = consistent_predict_index(predict_index)
        if (missing(realized.y)) {
            stop("You must provide a realization of the dependent variable in realized.y")
        }
        return(dens_yf(consistent.yf(realized.y, predict_index), 
            predict_index))
    }
    reslist$lps = function(realized.y, predict_index = NULL) {
        predict_index = consistent_predict_index(predict_index)
        if (missing(realized.y)) {
            stop("You must provide a realization of the dependent variable in realized.y")
        }
        yf = consistent.yf(realized.y, predict_index)
        if (ncol(yf) != 1) 
            stop("realized.y must have only one column!")
        yf.dens = dens_yf(yf, predict_index)
        return(-sum(log(yf.dens[!is.na(yf.dens)]))/length(yf))
    }
    reslist$plot = function(predict_index = NULL, addons = "eslz", 
        realized.y = NULL, addons.lwd = 1.5, ...) {
        dotargs = match.call(expand.dots = FALSE)$...
        xf_series = seq_len(rnew)
        predict_index = consistent_predict_index(predict_index)
        if (is.null(predict_index)) {
            predict_index = xf_series
        }
        else if (!all(predict_index %in% xf_series)) 
            stop(paste("predict_index needs to be an integer between 1 and ", 
                rnew, "!", sep = ""))
        if (!(is.null(realized.y))) {
            if (length(realized.y) != length(predict_index)) {
                stop("realized.y must be a vector with the same number of elements as rows in newdata (or predict_index)!")
            }
        }
        if (!is.null(realized.y)) 
            realized.y <- consistent.yf(realized.y, predict_index)
        calc_alldens()
        oldask = par()$ask
        plotnb = 0
        for (xf_index in predict_index) {
            doplot = !dlist[[xf_index]]$has.na
            plotnb = plotnb + doplot
            if (plotnb == 2) 
                par(ask = TRUE)
            dotargs = .adjustdots(dotargs, main = NULL, col = "steelblue4", 
                xlab = "Response variable")
            if (doplot) {
                eval(as.call(c(list(as.name("plot.preddens"), 
                  as.name("xf_index"), addons = as.name("addons"), 
                  yf.addons = as.name("realized.y"), addons.lwd = as.name("addons.lwd")), 
                  as.list(dotargs))))
            }
        }
        par(ask = oldask)
    }
    reslist$n = n
    reslist$nmodel = nmodel
    reslist$call = sys.call(0)
    class(reslist) = "pred.density"
    rm(betahat_elem1, bools, emptydens, firstbetas, firstbetas2, 
        linvres)
    return(reslist)
}
