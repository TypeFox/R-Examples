.adjustdots <-
function (dotargs, ...) 
{
    defargs = list(...)
    defargnames = names(defargs)
    dotargs = as.list(dotargs)
    if (is.null(dotargs)) {
        dotargs = list()
    }
    for (di in seq_len(length(defargs))) {
        if (!is.element(defargnames[[di]], names(dotargs))) {
            dotargs[[defargnames[[di]]]] <- defargs[[di]]
        }
    }
    return(dotargs)
}
.choose.gprior <-
function (g, N, K, return.g.stats = FALSE, yty = N, ...) 
{
    if (is.list(g)) {
        if (!all(c("gtype", "is.constant", "return.g.stats", 
            "lprobcalc") %in% names(g))) 
            stop("The provided g-prior list (in argument 'g') does not conform to the standards of a g-prior object.")
        if (!("g" %in% names(g))) {
            g$is.constant = FALSE
            g$g = NA
        }
        if (!("shrinkage.moments" %in% names(g))) {
            g$shrinkage.moments = ifelse((g$is.constant && is.numeric(g)), 
                g/(1 + g), 0)
        }
        if (!all(sapply(g$lprobcalc, is.function))) 
            stop("The slot 'lprobcalc' in the provided g-prior list (in argument 'g') does not conform to the standards of a g-prior object.")
        return(g)
    }
    if (is.function(g)) {
        return(g(g = g, return.g.stats = return.g.stats, N = N, 
            K = K, yty = yty, ...))
    }
    if (is.numeric(g)) {
        return(.gprior.constg.init(g = g, return.g.stats = return.g.stats, 
            N = N, K = K, yty = yty))
    }
    if (any(grep("EBL", g, ignore.case = TRUE))) {
        return(.gprior.eblocal.init(g = g, return.g.stats = return.g.stats, 
            N = N, K = K, yty = yty))
    }
    if (any(grep("hyper", g, ignore.case = TRUE))) {
        return(.gprior.hyperg.init(g = g, return.g.stats = return.g.stats, 
            N = N, K = K, yty = yty))
    }
    return(.gprior.constg.init(g = g, return.g.stats = return.g.stats, 
        N = N, K = K, yty = yty))
}
.choose.mprior <-
function (mpmode, mpparam, K, ..., fixed.pos = numeric(0)) 
{
    origargs = list(mpmode = mpmode, mpparam = mpparam)
    fixed.pos = {
        1:K
    }[fixed.pos]
    fixed.exist = as.logical(length(fixed.pos))
    fixk = length(fixed.pos)
    if (!(is.character(mpmode) || is.function(mpmode) || is.list(mpmode))) 
        stop("'mprior' parameter must be character! (or function/list)")
    if (is.function(mpmode) || is.list(mpmode)) {
        mpinfo = mpmode
    }
    else if (any(grep("fix", mpmode, ignore.case = TRUE))) {
        mpinfo = .mprior.fixedt.init
        if (is.numeric(mpparam)) 
            mpparam = mpparam[[1]] - fixk
    }
    else if (any(grep("unif", mpmode, ignore.case = TRUE))) {
        mpinfo = .mprior.uniform.init
    }
    else if (any(grep("custom", mpmode, ignore.case = TRUE))) {
        mpinfo = .mprior.customk.init
        if (fixed.exist && is.numeric(mpparam)) 
            if (length(mpparam) == K + 1) 
                mpparam = mpparam[(fixk + 1):length(mpparam)]
    }
    else if (any(grep("pip", mpmode, ignore.case = TRUE))) {
        mpinfo = .mprior.pip.init
        if (fixed.exist && is.numeric(mpparam)) 
            if (length(mpparam) == K) 
                mpparam = mpparam[-fixed.pos]
    }
    else {
        mpinfo = .mprior.randomt.init
        if (is.numeric(mpparam)) 
            mpparam = mpparam[[1]] - fixk
    }
    if (is.function(mpinfo)) 
        mpinfo = .fixedset.mprior(mpinfo, fullK = K, fixed.pos = fixed.pos, 
            K = NA, mpparam = mpparam, mpmode = mpmode, ...)
    if (!all(c("mp.mode", "mp.msize", "pmp") %in% names(mpinfo))) 
        stop("The provided custom-built model prior is deficient.")
    if (!("origargs" %in% names(mpinfo))) 
        mpinfo$origargs = origargs
    if (length(fixed.pos) > 0) 
        mpinfo$fixed.pos = fixed.pos
    class(mpinfo) <- c("mprior", class(mpinfo))
    return(mpinfo)
}
.constr.intmat <-
function (X, K) 
{
    intix = grep("#", colnames(X), fixed = TRUE)
    mPlus = diag(K)
    colnames(mPlus) <- colnames(X)
    for (jj in 1:length(intix)) {
        cix = intix[jj]
        mPlus[cix, unlist(strsplit(colnames(mPlus)[cix], "#", 
            fixed = TRUE))] = 1
    }
    return(mPlus)
}
.construct.arglist <-
function (funobj, envir = NULL) 
{
    namedlist = formals(funobj)
    argnames = names(namedlist)
    if (!is.environment(envir)) 
        envir = sys.frame(-1)
    for (argn in 1:length(namedlist)) {
        testval = as.logical(try(exists(argnames[argn], envir = envir), 
            silent = TRUE))
        if (is.na(testval)) 
            testval = FALSE
        if (testval) {
            namedlist[[argn]] = try(get(argnames[argn], envir = envir))
        }
    }
    return(namedlist)
}
.cor.topmod <-
function (tmo) 
{
    if (is.bma(tmo)) 
        tmo = tmo$topmod
    pmp.10 = pmp.bma(tmo, oldstyle = TRUE)
    if (nrow(pmp.10) == 1 | suppressWarnings(length(grep("error", 
        class(try(cor(pmp.10[, 1], pmp.10[, 2]), silent = TRUE)))))) {
        corr.pmp = NA
    }
    else {
        if (var(pmp.10[, 2]) == 0) 
            corr.pmp = NA
        else corr.pmp = cor(pmp.10[, 1], pmp.10[, 2])
    }
    return(corr.pmp)
}
.enum_fromindex <-
function (lindex) 
{
    lindex = lindex[[1]]
    if (lindex == 0) 
        return(FALSE)
    log2 = ceiling(log(lindex + 1, 2))
    return(as.logical((lindex + 2^((log2 - 1):0))%/%(2^(log2:1))%%2))
}
.enum_startend <-
function (iter = NA, start.value = 0, K = 1, maxk = K, fixed.pos = numeric(0)) 
{
    fixed.pos = {
        1:K
    }[fixed.pos]
    effk = K - length(fixed.pos)
    flexpos = {
        1:K
    }
    if (length(fixed.pos) > 0) 
        flexpos = {
            1:K
        }[-fixed.pos]
    start.value2 = 0
    if (length(start.value) == 1) {
        start.value2 = suppressWarnings(as.integer(start.value))
        if (any(is.na(start.value2)) | start.value2[[1]] < 0 | 
            start.value2[[1]] >= (2^effk - 1)) {
            start.value = 0
            start.value2 = 0
        }
    }
    else {
        start.value = 0
    }
    if (length(start.value) == 1) {
        start.value_cut = .enum_fromindex(start.value2)
        start.value = rep(1, K)
        start.value[flexpos] = c(numeric(effk - length(start.value_cut)), 
            start.value_cut)
    }
    if (K > maxk) {
        lastindex = 2^effk - 1 - sum(choose(effk, (maxk + 1):effk))
    }
    else lastindex = 2^effk - 1
    if (is.na(iter)) {
        iter = lastindex - start.value2
    }
    iter = min(iter, 2^effk - 1 - start.value2)
    return(list(start.value = start.value, iter = iter))
}
.f21_4hyperg <-
function (N, K, f21a, ltermbounds = c(200, 600, 1400, 3000)) 
{
    create.lterms = function(cc) lapply(as.list(ltermbounds), 
        function(x) (a + 0:x)/(cc + 0:x))
    a = (N - 1)/2
    cv = (f21a + 0:K)/2
    lterms = (lapply(cv, create.lterms))
    ltermbounds = c(0, ltermbounds[-length(ltermbounds)])
    return(list(calcit = function(z, k) {
        nbterms = sum(ceiling(abs({
            a - cv[k + 1]
        }/{
            1 - z
        }) * 1.5) >= ltermbounds)
        return(sum(cumprod(z * lterms[[k + 1]][[nbterms]])) + 
            1)
    }))
}
.f21simple <-
function (a, c, z) 
{
    f21o = .f21_4hyperg(2 * a + 1, 0, c * 2)
    f21o$calcit(z, 0)
}
.fixedset.mprior <-
function (mprior.function, fullK, fixed.pos = numeric(0), K = NA, 
    ...) 
{
    if (length(fixed.pos) == 0) 
        return(mprior.function(K = fullK, ...))
    fixed.pos = {
        1:fullK
    }[fixed.pos]
    flexpos = {
        1:fullK
    }[-fixed.pos]
    flexk = length(flexpos)
    mprior = mprior.function(K = flexk, ...)
    fixk = length(fixed.pos)
    mpl = list(mp.mode = mprior$mp.mode, mp.msize = mprior$mp.msize + 
        fixk, pmp = function(ki, mdraw, ...) {
        return(mprior$pmp(ki = ki - fixk, mdraw = mdraw[flexpos]))
    }, mp.Kdist = c(numeric(fixk), mprior$mp.Kdist))
    return(mpl)
}
.fixedset.sampler <-
function (sampler.function, fullK, fixed.pos = numeric(0), ...) 
{
    if (length(fixed.pos) == 0) 
        return(sampler.function)
    fixed.pos = {
        1:fullK
    }[fixed.pos]
    flexpos = {
        1:fullK
    }[-fixed.pos]
    flexk = length(flexpos)
    outdraw = rep(1, fullK)
    outfun = function(molddraw = molddraw, K = flexk, ...) {
        flexdraw = sampler.function(molddraw = molddraw[flexpos], 
            K = flexk, ...)
        outdraw[flexpos] = flexdraw[["mnewdraw"]]
        addi = flexdraw[["addi"]]
        dropi = flexdraw[["dropi"]]
        if (is.numeric(addi) || is.numeric(dropi)) {
            if (addi > 0) 
                addi = flexpos[addi]
            else addi = 0
            if (dropi > 0) 
                dropi = flexpos[dropi]
            else dropi = 0
        }
        return(list(mnewdraw = outdraw, positionnew = {
            1:fullK
        }[as.logical(outdraw)], addi = addi, dropi = dropi))
    }
    return(outfun)
}
.fls.samp <-
function (molddraw = molddraw, K = K, ..., maxk = Inf, oldk = 0) 
{
    indch <- ceiling(runif(1, 0, K))
    bdropit <- as.logical(molddraw[[indch]])
    if (oldk == maxk) 
        if (!bdropit) {
            indch = (1:K)[molddraw == 1][[ceiling(runif(1, 0, 
                sum(molddraw)))]]
            bdropit = molddraw[[indch]]
        }
    if (bdropit) {
        addvar <- 0
        dropvar <- indch
        molddraw[[indch]] <- 0
    }
    else {
        addvar <- indch
        dropvar <- 0
        molddraw[[indch]] <- 1
    }
    positionnew <- {
        1:K
    }[molddraw == 1]
    return(list(mnewdraw = molddraw, positionnew = positionnew, 
        addi = addvar, dropi = dropvar))
}
.fls.samp.int <-
function (molddraw = molddraw, K = K, mPlus = mPlus, maxk = Inf, 
    oldk = 0) 
{
    indch = ceiling(runif(1, 0, 1) * K)
    if (molddraw[indch] == 1) {
        mnewdraw = as.numeric(molddraw > mPlus[, indch])
        dropvar = (1:K)[xor(molddraw, mnewdraw)]
        addvar = 0
    }
    else {
        mnewdraw = as.numeric(molddraw | mPlus[indch, ])
        addvar = (1:K)[xor(molddraw, mnewdraw)]
        dropvar = 0
    }
    positionnew = which(mnewdraw == 1)
    if (length(positionnew) > maxk) {
        return(.fls.samp.int(molddraw = molddraw, K = K, mPlus = mPlus, 
            maxk, oldk))
    }
    else {
        return(list(mnewdraw = mnewdraw, positionnew = positionnew, 
            addi = addvar, dropi = dropvar))
    }
}
.getpolycoefs <-
function (polyroots) 
{
    if (length(polyroots) == 1) 
        return(c(1, polyroots))
    restterms = .getpolycoefs(polyroots[-1])
    c(restterms, 0) + c(0, polyroots[1] * restterms)
}
.gprior.constg.init <-
function (g = NA, return.g.stats = TRUE, N = N, K = K, yty = 1, 
    null.lik = NA, ...) 
{
    gg = NULL
    if (!(is.character(g) || is.numeric(g))) 
        g = "UIP"
    if (any(grep("BRIC", g, ignore.case = TRUE)) | any(grep("FLS", 
        g, ignore.case = TRUE))) {
        if (N <= (K^2)) {
            gg = (K^2)
        }
        else {
            gg = N
        }
        gtype = "BRIC"
    }
    if (any(grep("RIC", g, ignore.case = TRUE)) && (!any(grep("BRIC", 
        g, ignore.case = TRUE)))) {
        gg = (K^2)
        gtype = "RIC"
    }
    if (any(grep("HQ", g, ignore.case = TRUE)) | any(grep("Hannan", 
        g, ignore.case = TRUE))) {
        gg = (log(N))^3
        gtype = "Hannan-Quinn"
    }
    if (is.numeric(g)) {
        gg = g
        gtype = "numeric"
    }
    if (is.null(gg)) {
        if (!(any(grep("UIP", g, ignore.case = TRUE)) | any(grep("BIC", 
            g, ignore.case = TRUE)))) 
            warning("The provided g prior could not be identified. Therefore the default g prior (UIP) has been selected.")
        gg = N
        gtype = "UIP"
    }
    gprior.info = list(gtype = gtype, is.constant = TRUE, return.g.stats = return.g.stats, 
        shrinkage.moments = gg/(gg + 1), g = gg)
    g = gg
    if (!is.numeric(null.lik)) {
        null.lik = {
            1 - N
        }/2 * log(yty)
    }
    g2 = g/{
        g + 1
    }
    l1g = log(1 + g)
    g2sq = g2^2
    n1 = N - 1
    gprior.info$lprobcalc <- list(just.loglik = function(ymy, 
        k, ...) {
        return(0.5 * {
            {
                n1 - k
            } * l1g - n1 * log(g * ymy + yty)
        })
    }, lprob.all = function(ymy, k, bhat, diag.inverse, ...) {
        b1new = g2 * bhat
        return(list(lprob = 0.5 * {
            {
                n1 - k
            } * l1g - n1 * log(g * ymy + yty)
        }, b1new = b1new, b2new = {
            {
                yty/g + ymy
            } * g2sq/{
                N - 3
            }
        } * diag.inverse + b1new^2, otherstats = numeric(0)))
    })
    class(gprior.info) <- c("gprior", class(gprior.info))
    return(gprior.info)
}
.gprior.eblocal.init <-
function (g = NA, return.g.stats = TRUE, N = N, K = K, yty = 1, 
    null.lik = NA, ...) 
{
    gprior.info = list(gtype = "EBL", is.constant = FALSE, return.g.stats = return.g.stats, 
        shrinkage.moments = numeric(1), g = NA)
    if (!is.numeric(null.lik)) {
        null.lik = (1 - N)/2 * log(yty)
    }
    ymy.current = -1
    k.current = -1
    loglik = null.lik
    Fstat = numeric(0)
    g2 = 0
    if (return.g.stats) {
        otherstats = numeric(1)
    }
    else {
        otherstats = numeric(0)
    }
    gprior.info$lprobcalc <- list(just.loglik = function(ymy, 
        k, ...) {
        ymy.current <<- ymy
        k.current <<- k
        if (k == 0) {
            return(null.lik)
        }
        Fstat <<- (N - k - 1)/k * (yty - ymy)/ymy
        if (Fstat > 1) {
            g0 <- 1/(Fstat - 1)
            g2 <<- 1/(g0 + 1)
        } else {
            g2 <<- 0
            return(null.lik)
        }
        lFstat = log(Fstat)
        lg02 = -lFstat
        loglik <<- 0.5 * {
            k * lg02 - {
                N - 1
            } * {
                log(ymy + g0 * yty) + {
                  log(Fstat - 1) - lFstat
                }
            }
        }
        return(loglik)
    }, lprob.all = function(ymy, k, bhat, diag.inverse, ...) {
        if (k == 0) {
            return(list(lprob = null.lik, b1new = numeric(0), 
                b2new = numeric(0), otherstats = otherstats))
        }
        if ((ymy != ymy.current) | (k != k.current)) {
            Fstat <<- {
                N - k - 1
            }/k * (yty - ymy)/ymy
            if (Fstat > 1) {
                g0 = 1/{
                  Fstat - 1
                }
                g2 <<- 1/(g0 + 1)
                lFstat = log(Fstat)
                lg02 = -lFstat
                loglik <<- 0.5 * {
                  k * lg02 - {
                    N - 1
                  } * {
                    log(ymy + g0 * yty) + {
                      log(Fstat - 1) - lFstat
                    }
                  }
                }
            } else {
                g0 = 0
                g2 <<- 0
                loglik <<- null.lik
            }
        }
        if (return.g.stats) {
            otherstats = g2
        }
        b1new = g2 * bhat
        if (g2 > 0) {
            b2new = {
                {
                  (1/g2 - 1) * yty + ymy
                } * {
                  g2^2
                }/{
                  N - 3
                }
            } * diag.inverse + b1new^2
        } else {
            b2new = numeric(k)
        }
        return(list(lprob = loglik, b1new = b1new, b2new = b2new, 
            otherstats = otherstats))
    })
    class(gprior.info) <- c("gprior", class(gprior.info))
    return(gprior.info)
}
.gprior.hyperg.init <-
function (g = NA, return.g.stats = TRUE, N = N, K = K, yty = 1, 
    null.lik = NA, ...) 
{
    if (!is.character(g)) 
        g = "hyper"
    if (any(grep("=", g))) {
        f21a = suppressWarnings(as.numeric(unlist(strsplit(g, 
            "="))[2]))
        if (!is.numeric(f21a) | is.na(f21a)) {
            f21a.char = suppressWarnings(as.character(unlist(strsplit(g, 
                "="))[2]))
            if (any(grep("bric", f21a.char, ignore.case = TRUE))) {
                f21a = 2 + 2/max(N, K^2)
            }
            else if (any(grep("uip", f21a.char, ignore.case = TRUE))) {
                f21a = 2 + 2/N
            }
            else {
                warning("You did not supply a proper 'a' parameter for the hyper g prior (like e.g. the format g='hyperg=3.1' or g='hyper=UIP') - thus set to default value 'hyper=UIP' instead.")
                f21a = 2 + 2/N
            }
        }
        else {
            if (f21a <= 2 | f21a > 4) {
                f21a = 2 + 2/N
                warning("You provided an 'a' parameter for the hyper g prior that is not element of (2,4]. I chose the default value 'hyper=UIP' instead.")
            }
        }
    }
    else {
        f21a = 2 + 2/N
    }
    gprior.info = list(gtype = "hyper", is.constant = FALSE, 
        return.g.stats = return.g.stats, shrinkage.moments = numeric(2), 
        g = NA, hyper.parameter = f21a)
    if (!is.numeric(null.lik)) {
        null.lik = {
            1 - N
        }/2 * log(yty)
    }
    gmoments = numeric(2)
    N12 = {
        N - 1
    }/2
    la2 = log(f21a - 2)
    log.lik = null.lik
    ymy.current = -1
    k.current = -1
    intconstinv = f21a - 2
    f21o = .f21_4hyperg(N, K, f21a)
    gprior.info$lprobcalc <- list(just.loglik = function(ymy, 
        k, ...) {
        if (k == 0) {
            return(null.lik)
        }
        ymy.current <<- ymy
        k.current <<- k
        intconstinv <<- {
            k + f21a - 2
        }/f21o[["calcit"]](1 - ymy/yty, k)
        if (intconstinv < 0) {
            intconstinv <<- k + f21a - 2
        }
        log.lik <<- null.lik + la2 - log(intconstinv)
        return(log.lik)
    }, lprob.all = function(ymy, k, bhat, diag.inverse, ...) {
        if (k == 0) {
            return(list(lprob = null.lik, b1new = numeric(0), 
                b2new = numeric(0), otherstats = c(2/f21a, 8/f21a/(f21a + 
                  2))))
        }
        N3 = N - 3
        ka2 = k + f21a - 2
        R2 = 1 - ymy/yty
        if ((ymy != ymy.current) | (k != k.current)) {
            intconstinv <<- ka2/f21o[["calcit"]](R2, k)
            log.lik <<- null.lik + la2 - log(intconstinv)
        }
        g2hyper = {
            intconstinv - ka2 + N3 * R2
        }/{
            R2 * {
                N3 - ka2
            }
        }
        gbetavar = {
            {
                1 + 2/N3 * R2/{
                  1 - R2
                }
            } * intconstinv + {
                N3 - 2
            } * R2 - ka2
        } * N3 * {
            1 - R2
        }/{
            N3 - ka2
        }/{
            N3 - ka2 - 2
        }/R2 * yty/N3 * diag.inverse
        if (return.g.stats) {
            ka = ka2 + 2
            Eg22 = {
                {
                  {
                    N3 - 2
                  } * R2 - ka
                } * intconstinv + {
                  N3 * R2 - ka2
                }^2 - 2 * {
                  N3 * R2^2 - ka2
                }
            }/R2^2/{
                N3 - ka2
            }/{
                N3 - ka
            }
            gmoments = c(g2hyper, Eg22)
        }
        return(list(lprob = log.lik, b1new = g2hyper * bhat, 
            b2new = gbetavar + g2hyper^2 * bhat^2, otherstats = gmoments))
    })
    class(gprior.info) <- c("gprior", class(gprior.info))
    return(gprior.info)
}
.hexcode.binvec.convert <-
function (length.of.binvec) 
{
    if (length(length.of.binvec) > 1) 
        length.of.binvec = length(length.of.binvec)
    addpositions = 4 - length.of.binvec%%4
    positionsby4 = (length.of.binvec + addpositions)/4
    hexvec = c(0:9, "a", "b", "c", "d", "e", "f")
    hexcodelist = list(`0` = numeric(4), `1` = c(0, 0, 0, 1), 
        `2` = c(0, 0, 1, 0), `3` = c(0, 0, 1, 1), `4` = c(0, 
            1, 0, 0), `5` = c(0, 1, 0, 1), `6` = c(0, 1, 1, 0), 
        `7` = c(0, 1, 1, 1), `8` = c(1, 0, 0, 0), `9` = c(1, 
            0, 0, 1), a = c(1, 0, 1, 0), b = c(1, 0, 1, 1), c = c(1, 
            1, 0, 0), d = c(1, 1, 0, 1), e = c(1, 1, 1, 0), f = c(1, 
            1, 1, 1))
    return(list(as.hexcode = function(binvec) {
        incl = c(numeric(addpositions), binvec)
        dim(incl) = c(4, positionsby4)
        return(paste(hexvec[crossprod(incl, 2L^(3:0)) + 1], collapse = ""))
    }, as.binvec = function(hexcode) {
        return(unlist(hexcodelist[unlist(strsplit(hexcode, "", 
            fixed = TRUE), recursive = FALSE, use.names = FALSE)], 
            recursive = FALSE, use.names = FALSE)[-(1:addpositions)])
    }))
}
.index.bma <-
function (x, i, ...) 
{
    x$topmod <- x$topmod[i]
    return(x)
}
.index.topmod <-
function (x, i, ...) 
{
    tm = x
    idx = i
    if (any(is.na(suppressWarnings(as.integer(idx))))) 
        idx = 1:length(tm$lik())
    if (length(tm$betas_raw()) > 1) {
        bbeta = TRUE
        bet = as.vector(tm$betas()[, idx])
        bet = bet[bet != 0]
    }
    else {
        bbeta = FALSE
        bet = numeric(0)
    }
    if (length(tm$betas2_raw()) > 1) {
        bbeta2 = TRUE
        bet2 = as.vector(tm$betas2()[, idx])
        bet2 = bet2[bet2 != 0]
    }
    else {
        bbeta2 = FALSE
        bet2 = numeric(0)
    }
    fixvec = tm$fixed_vector()
    if (!length(as.vector(fixvec))) 
        fixvec = numeric(0)
    else fixvec = as.vector(t(fixvec[, idx]))
    .top10(nmaxregressors = tm$nregs, nbmodels = tm$nbmodels, 
        bbeta = bbeta, lengthfixedvec = nrow(tm$fixed_vector()), 
        bbeta2 = bbeta2, inivec_lik = tm$lik()[idx], inivec_bool = tm$bool()[idx], 
        inivec_count = tm$ncount()[idx], inivec_vbeta = bet, 
        inivec_vbeta2 = bet2, inivec_veck = tm$kvec_raw()[idx], 
        inivec_fixvec = fixvec)
}
.iterenum <-
function (molddraw = numeric(0), K = length(molddraw), ...) 
{
    even.lead1 = {
        1:K
    }[!{
        cumsum(molddraw)%%2
    }]
    i = even.lead1[length(even.lead1)]
    molddraw[i] = !molddraw[i]
    addi = molddraw[i] * i
    dropi = {
        !molddraw[i]
    } * i
    return(list(mnewdraw = molddraw, positionnew = {
        1:K
    }[as.logical(molddraw)], addi = addi, dropi = dropi))
}
.iterenum.bone <-
function (molddraw = numeric(0), maxk = Inf) 
{
    even.lead1 = ((1:length(molddraw))[!(cumsum(molddraw)%%2)])
    i = even.lead1[length(even.lead1)]
    molddraw[i] = !molddraw[i]
    if (sum(molddraw) > maxk) 
        return(.iterenum.bone(molddraw, maxk))
    else return(molddraw)
}
.iterenum.KgtN <-
function (molddraw = numeric(0), maxk = Inf, oldk = 0, ...) 
{
    mnewdraw = .iterenum.bone(molddraw = molddraw, maxk)
    addi = (1:length(mnewdraw))[molddraw < mnewdraw]
    if (length(addi) == 0) 
        addi = 0
    dropi = (1:length(mnewdraw))[molddraw > mnewdraw]
    if (length(dropi) == 0) 
        dropi = 0
    return(list(mnewdraw = mnewdraw, positionnew = (1:length(mnewdraw))[as.logical(mnewdraw)], 
        addi = addi, dropi = dropi))
}
.lprob.constg.init <-
function (...) 
{
    gpo = .gprior.constg.init(...)
    return(gpo$lprobcalc)
}
.lprob.eblocal.init <-
function (...) 
{
    gpo = .gprior.eblocal.init(...)
    return(gpo$lprobcalc)
}
.lprob.hyperg.init <-
function (...) 
{
    gpo = .gprior.hyperg.init(...)
    return(gpo$lprobcalc)
}
.mprior.customk.init <-
function (K, mpparam, ...) 
{
    if (any(is.na(mpparam))) 
        mpparam = rep(0.5, K)
    if (!is.numeric(mpparam)) 
        stop("For custom model size priors, you need to provide a K+1 vector with positive elements for argument 'mprior.size'.")
    mpparam = as.vector(mpparam)
    if (!((length(mpparam) == (K + 1)) & all(mpparam > 0))) {
        stop("For custom model size priors, you need to provide a K+1 vector with positive elements for argument 'mprior.size'.")
    }
    mpkvec = log(mpparam)
    return(list(mp.mode = "custom", mp.msize = sum(choose(K, 
        0:K) * mpparam * {
        0:K
    })/sum(choose(K, 0:K) * mpparam), pmp = function(ki, ...) {
        return(mpkvec[[ki + 1]])
    }, mp.Kdist = choose(K, 0:K) * mpparam/sum(choose(K, 0:K) * 
        mpparam)))
}
.mprior.fixedt.init <-
function (K, mpparam, ...) 
{
    if (is.na(mpparam[1])) 
        mpparam <- K/2
    if ((mpparam[[1]] >= K) & (length(mpparam) == 1)) {
        warning("Submitted prior model size is >= than the nr. of   regressors\n, used K/2 instead\n\n")
        mpparam <- K/2
    }
    m = mpparam[[1]]
    return(list(mp.mode = "fixed", mp.msize = m, pmp = function(ki, 
        ...) {
        post.odds1 = ki * log(m/K) + {
            K - ki
        } * log(1 - m/K)
        return(post.odds1)
    }, mp.Kdist = dbinom(x = 0:K, size = K, prob = m/K, log = FALSE)))
}
.mprior.pip.init <-
function (K, mpparam, ...) 
{
    if (any(is.na(mpparam))) 
        mpparam = rep(0.5, K)
    if (!is.numeric(mpparam)) 
        stop("For prior inclusion probabilites, you need to provide a K vector with elements between 0 and 1 for argument 'mprior.size'.")
    mpparam = as.vector(mpparam)
    if (!((length(mpparam) == K) & all(mpparam > 0) & all(mpparam <= 
        1))) 
        stop("For prior inclusion probabilites, you need to provide a K vector with elements between 0 and 1 for argument 'mprior.size'.")
    if (any(mpparam == 1L)) 
        warning("Prior Inclsuion Prob. = 1 are impractical. Try using the argument fixed.reg")
    inclfacts = log(mpparam/(1 - mpparam))
    return(list(mp.mode = "pip", mp.msize = sum(mpparam), pmp = function(mdraw, 
        ...) {
        return(sum(inclfacts[as.logical(mdraw)]))
    }, mp.Kdist = .getpolycoefs(mpparam/{
        1 - mpparam
    }) * prod(1 - mpparam)))
}
.mprior.randomt.init <-
function (K, mpparam, ...) 
{
    if (is.na(mpparam[1])) 
        mpparam <- K/2
    if ((mpparam[[1]] >= K) & (length(mpparam) == 1)) {
        warning("Submitted prior model size is >= than the nr. of   regressors\n, used K/2 instead\n\n")
        mpparam <- K/2
    }
    m = mpparam[[1]]
    vecofpriors = lgamma(1 + 0:K) + lgamma({
        K - m
    }/m + K - 0:K)
    beta.bin = function(a = 1, b = (K - m)/m, K = K, w = 0:K) {
        lgamma(a + b) - {
            lgamma(a) + lgamma(b) + lgamma(a + b + K)
        } + log(choose(K, w)) + lgamma(a + w) + lgamma(b + K - 
            w)
    }
    return(list(mp.mode = "random", mp.msize = m, pmp = function(ki, 
        ...) {
        return(vecofpriors[[ki + 1]])
    }, mp.Kdist = exp(beta.bin(a = 1, b = {
        K - m
    }/m, K = K, w = 0:K))))
}
.mprior.uniform.init <-
function (K, ...) 
{
    return(list(mp.mode = "uniform", mp.msize = K/2, pmp = function(...) return(0), 
        mp.Kdist = exp(lchoose(K, 0:K) - K * log(2))))
}
.ols.terms2 <-
function (positions, yty, k = NULL, N = N, K = K, XtX.big = XtX.big, 
    Xty.big = Xty.big, ...) 
{
    syminv <- function(symmat, ndim = ncol(symmat)) {
        if (!is.matrix(symmat)) {
            symmat = as.matrix(symmat)
        }
        return(chol2inv(chol(symmat), size = ndim))
    }
    if (is.null(k)) 
        k = length(positions)
    XtXinv.return = numeric(0)
    if (sum(k) == 0) {
        Xty = numeric(0)
        XtXinv = matrix(0, 0, 0)
        bhat = numeric(0)
        ymy = yty
        positions = 0
    }
    else {
        XtX <- XtX.big[positions, positions, drop = FALSE]
        Xty <- Xty.big[positions]
        XtXinv <- syminv(XtX, ndim = k)
        bhat <- crossprod(XtXinv, Xty)
        ymy <- yty - crossprod(Xty, bhat)[[1]]
    }
    return(list(full.results = function() {
        return(list(ymy = ymy, bhat = bhat, diag.inverse = XtXinv[1:k + 
            0:(k - 1) * k]))
    }, child.ymy = function(addix = 0, dropix = 0, ...) {
        if (!any(as.logical(c(addix, dropix)))) {
            return(ymy)
        }
        if (all(as.logical(c(addix, dropix)))) {
            jhere = {
                1:k
            }[positions == dropix]
            poshere = positions[-jhere]
            Xj = XtXinv[, jhere]
            Xtxi = XtX.big[poshere, addix]
            bxlessj = crossprod(XtXinv, XtX.big[positions, addix]) - 
                Xj * XtX.big[addix, dropix]
            bhatx = bxlessj[-jhere] - Xj[-jhere] * bxlessj[jhere]/Xj[jhere]
            child.ymy = ymy + bhat[jhere]^2/Xj[jhere] - {
                Xty.big[addix] - crossprod(Xty.big[poshere], 
                  bhatx)[[1]]
            }^2/{
                XtX.big[addix, addix] - crossprod(bhatx, Xtxi)[[1]]
            }
            return(child.ymy)
        } else {
            if (addix == 0) {
                jhere = {
                  1:k
                }[positions == dropix]
                child.ymy = ymy + bhat[jhere]^2/XtXinv[jhere, 
                  jhere]
                return(child.ymy)
            } else {
                Xtxi = XtX.big[positions, addix]
                bhatx = crossprod(XtXinv, Xtxi)[, 1]
                child.ymy = ymy - {
                  Xty.big[addix] - crossprod(bhatx, Xty)[[1]]
                }^2/{
                  XtX.big[addix, addix] - crossprod(bhatx, Xtxi)[[1]]
                }
                return(child.ymy)
            }
        }
    }, mutate = function(addix = 0, dropix = 0, newpos = numeric(0), 
        newk = 0, ...) {
        if (newk == 0) {
            XtXinv <<- matrix(0, 0, 0)
            Xty <<- numeric(0)
        } else {
            if (newk < 7 | addix[[1]] != 0 | length(c(dropix, 
                addix)) > 2) {
                Xty <<- Xty.big[newpos]
                XtXinv <<- syminv(XtX.big[newpos, newpos, drop = FALSE], 
                  ndim = newk)
            } else {
                if (dropix[1] > 0) {
                  jhere = sum(positions <= dropix)
                  Xty <<- Xty[-jhere]
                  Xj = XtXinv[, jhere]
                  XtXinv <<- {
                    XtXinv - tcrossprod(Xj/Xj[jhere], Xj)
                  }[-jhere, -jhere]
                } else {
                  jhere = sum(positions < addix) + 1
                  Xtxx = XtX.big[addix, newpos]
                  Xtx = Xtxx[-jhere]
                  Xty <<- Xty.big[newpos]
                  bhatx = crossprod(XtXinv, Xtx)[, 1]
                  bhatxadj = c(bhatx[0:(jhere - 1)], -1, bhatx[jhere:k])
                  if (jhere == newk) bhatxadj = bhatxadj[-(jhere + 
                    1:2)]
                  newinv = tcrossprod(bhatxadj, bhatxadj/(Xtxx[jhere] - 
                    crossprod(Xtx, bhatx)[[1]]))
                  newinv[-jhere, -jhere] = newinv[-jhere, -jhere] + 
                    XtXinv
                  XtXinv <<- newinv
                }
            }
        }
        positions <<- newpos
        k <<- newk
        bhat <<- crossprod(XtXinv, Xty)[, 1]
        ymy <<- yty - crossprod(Xty, bhat)[[1]]
        return(list(ymy = ymy, bhat = bhat, diag.inverse = XtXinv[1:k + 
            0:{
                k - 1
            } * k]))
    }, return.inverse = function() XtXinv, ymy = ymy, bhat = bhat, 
        diag.inverse = XtXinv[1:k + 0:{
            k - 1
        } * k]))
}
.post.beta.draws <-
function (topmods, reg.names, moment2 = FALSE) 
{
    if (moment2) 
        beta.draws = as.matrix(topmods$betas2())
    else beta.draws = as.matrix(topmods$betas())
    if (sum(beta.draws) == 0) {
        stop("The tompod object provided does not have saved betas. cf. bbeta argument in function topmod")
    }
    if (nrow(beta.draws) != length(reg.names)) {
        rownames(beta.draws) = c(reg.names, "W-Index")
    }
    else {
        rownames(beta.draws) = c(reg.names)
    }
    beta.names = topmods$bool()
    if (length(which(beta.names == "0")) > 0) {
        colnames(beta.draws) = beta.names[-c(which(beta.names == 
            "0"))]
    }
    else {
        colnames(beta.draws) = beta.names
    }
    return(beta.draws)
}
.post.calc <-
function (gprior.info, add.otherstats, k.vec, null.count, X.data, 
    topmods, b1mo, b2mo, iter, burn, inccount, models.visited, 
    K, N, msize, timed, cumsumweights = NA, mcmc = "bd", possign = NA) 
{
    postad.k.vec <- function(k.vec, null.count) c(null.count, 
        k.vec)
    postad.gprior.info <- function(gprior.info, add.otherstats = numeric(0), 
        cumsumweights = 1) {
        if (gprior.info$return.g.stats) {
            if (length(add.otherstats) > 0) {
                gprior.info$shrinkage.moments = add.otherstats/cumsumweights
            }
            else {
                gprior.info$shrinkage.moments = 1/(1 + 1/gprior.info$g)
            }
        }
        return(gprior.info)
    }
    postad.reg.names <- function(X.data) {
        if (is.null(colnames(X.data)[-1]) || colnames(X.data)[-1] == 
            "") {
            reg.names <- paste("beta", 1:K)
        }
        else {
            reg.names = colnames(X.data)[-1]
        }
        return(reg.names)
    }
    gprior.info = postad.gprior.info(gprior.info, add.otherstats, 
        cumsumweights)
    k.vec = postad.k.vec(k.vec, null.count)
    cons = .post.constant(X.data, b1mo/cumsumweights)
    pmp.10 = pmp.bma(topmods, oldstyle = TRUE)
    if (nrow(pmp.10) == 1 | suppressWarnings(length(grep("error", 
        class(try(cor(pmp.10[, 1], pmp.10[, 2]), silent = TRUE)))))) {
        corr.pmp = NA
    }
    else {
        if (var(pmp.10[, 2]) == 0) 
            corr.pmp = NA
        else corr.pmp = cor(pmp.10[, 1], pmp.10[, 2])
    }
    if (is.na(possign[[1]])) 
        possign = numeric(K)
    info.object = list(iter = iter, burn = burn, inccount = inccount, 
        models.visited = models.visited, b1mo = b1mo, b2mo = b2mo, 
        add.otherstats = add.otherstats, cumsumweights = cumsumweights, 
        K = K, N = N, corr.pmp = corr.pmp, msize = msize, timed = timed, 
        k.vec = k.vec, cons = cons, pos.sign = possign)
    reg.names = postad.reg.names(X.data)
    return(list(info = info.object, k.vec = k.vec, cons = cons, 
        gprior.info = gprior.info, pmp.10 = pmp.10, reg.names = reg.names))
}
.post.constant <-
function (X.data, Ebeta) 
{
    Xmeans = colMeans(X.data)
    cons = Xmeans[1] - crossprod(Ebeta, Xmeans[-1])
    return(as.vector(cons))
}
.post.estimates <-
function (b1mo = NULL, b2mo = NULL, cumsumweights = NULL, inccount = NULL, 
    topmods = NULL, X.data = NULL, reg.names = NULL, pos.sign = NULL, 
    exact = FALSE, order.by.pip = TRUE, include.constant = FALSE, 
    incl.possign = TRUE, std.coefs = FALSE, condi.coef = FALSE) 
{
    idx = 1:(length(b1mo))
    if (exact) {
        lt1 = topmods$lik() - max(topmods$lik())
        exact.pmp = as.vector(exp(lt1)/sum(exp(lt1)))
        pip = as.vector(topmods$bool_binary() %*% exact.pmp)
        idx = 1:(length(pip))
        betas = topmods$betas()
        betas2 = topmods$betas2()
        K = nrow(betas)
        Eb1 = tcrossprod(betas, t(exact.pmp))[, 1]
        Eb2 = tcrossprod(betas2, t(exact.pmp))[, 1]
        Ebsd = sqrt(Eb2 - Eb1^2)
        possign = round(tcrossprod((betas > 0), t(exact.pmp))[, 
            1]/pip, 8)
        possign[is.nan(possign)] = NA
    }
    else {
        pip = inccount/cumsumweights
        Eb1 = b1mo/cumsumweights
        Eb2 = b2mo/cumsumweights
        Ebsd = sqrt(Eb2 - Eb1^2)
        possign = round(pos.sign/inccount, 8)
        possign[is.nan(possign)] = NA
    }
    if (include.constant) 
        constterm = .post.constant(X.data, Eb1)
    if (condi.coef) {
        Eb1 = Eb1/pip
        Eb2 = Eb2/pip
        Ebsd = sqrt(Eb2 - Eb1^2)
        Eb1[is.nan(Eb1)] = 0
        Ebsd[is.nan(Ebsd)] = 0
    }
    if (std.coefs) {
        sddata = apply(as.matrix(X.data), 2, stats::sd)
        Eb1 = Eb1/sddata[1] * sddata[-1]
        Ebsd = Ebsd/sddata[1] * sddata[-1]
        if (include.constant) 
            constterm = constterm/sddata[1]
    }
    if (incl.possign) {
        post.mean <- cbind(pip, Eb1, Ebsd, possign, idx)
        rownames(post.mean) <- reg.names
        colnames(post.mean) <- c("PIP", "Post Mean", "Post SD", 
            "Cond.Pos.Sign", "Idx")
    }
    else {
        post.mean <- cbind(pip, Eb1, Ebsd, idx)
        rownames(post.mean) <- reg.names
        colnames(post.mean) <- c("PIP", "Post Mean", "Post SD", 
            "Idx")
    }
    if (order.by.pip) {
        post.mean <- post.mean[order(-post.mean[, 1]), ]
    }
    if (include.constant) {
        constrow = matrix(c(1, constterm, NA, rep(NA, incl.possign), 
            0), 1)
        rownames(constrow) = "(Intercept)"
        post.mean = rbind(post.mean, constrow)
    }
    return(post.mean)
}
.post.topmod.bma <-
function (topmods, reg.names = numeric(0)) 
{
    pmps = pmp.bma(topmods)
    if (is.bma(topmods)) {
        reg.names = topmods$reg.names
        topmods = topmods$topmod
    }
    rbind(.post.topmod.includes(topmods, reg.names), t(pmps))
}
.post.topmod.includes <-
function (topmods, reg.names) 
{
    topmod = topmods$bool_binary()
    colnames(topmod) <- topmods$bool()
    rownames(topmod) = reg.names
    return(topmod)
}
.quantile.density <-
function (x, probs = seq(0.25, 0.75, 0.25), names = TRUE, normalize = TRUE, 
    ...) 
{
    my.quantile.density = function(x, probs, names, normalize, 
        ...) {
        ycs = (cumsum(x$y) - (x$y - x$y[[1]])/2) * diff(x$x[1:2])
        if (normalize) 
            ycs = ycs/(ycs[[length(ycs)]])
        xin = x$x
        maxi = length(ycs)
        qqs = sapply(as.list(probs), function(qu) {
            iii = sum(ycs <= qu)
            if (iii == maxi) 
                return(Inf)
            else if (iii == 0L) 
                return(-Inf)
            else {
                return(xin[[iii + 1]] + ((ycs[[iii + 1]] - qu)/(ycs[[iii + 
                  1]] - ycs[[iii]])) * (xin[[iii]] - xin[[iii + 
                  1]]))
            }
        })
        if (as.logical(names)) 
            names(qqs) = paste(format(100 * probs, trim = TRUE, 
                digits = max(2L, getOption("digits"))), "%", 
                sep = "")
        return(qqs)
    }
    probs = as.vector(probs)
    if (is.element("density", class(x))) 
        return(my.quantile.density(x = x, probs = probs, names = names, 
            normalize = normalize))
    if (!all(sapply(x, function(dd) is.element("density", class(dd))))) 
        stop("x needs to be a density or list of densities")
    if (length(x) == 1L) 
        return(my.quantile.density(x = x[[1]], probs = probs, 
            names = names, normalize = normalize))
    qout = sapply(x, my.quantile.density, probs = probs, names = FALSE, 
        normalize = normalize)
    if (!is.matrix(qout)) {
        if (length(probs) > 1) 
            return(qout)
        qout = as.matrix(qout)
    }
    else qout = t(qout)
    if (as.logical(names)) 
        colnames(qout) = paste(format(100 * probs, trim = TRUE, 
            digits = max(2L, getOption("digits"))), "%", sep = "")
    return(qout)
}
.rev.jump <-
function (molddraw = molddraw, K = K, ..., maxk = Inf, oldk = 0) 
{
    rev.idx = ceiling(runif(1, 0, 2))
    if (rev.idx == 1) {
        birth.death = .fls.samp(molddraw = molddraw, K = K, maxk = maxk, 
            oldk = oldk)
        mnewdraw = birth.death[["mnewdraw"]]
        positionnew = birth.death[["positionnew"]]
        addvar = birth.death[["addi"]]
        dropvar = birth.death[["dropi"]]
    }
    if (rev.idx == 2) {
        var.in = (1:K)[as.logical(molddraw)]
        var.out = (1:K)[!as.logical(molddraw)]
        var.in.rand = ceiling(length(var.in) * runif(1, 0, 1))
        addvar = var.out[ceiling(length(var.out) * runif(1, 0, 
            1))]
        dropvar = var.in[var.in.rand]
        mnewdraw = molddraw
        mnewdraw[addvar] = 1
        mnewdraw[dropvar] = 0
        positionnew = (1:K)[as.logical(mnewdraw)]
        dropvar = max(dropvar, 0)
        addvar = max(addvar, 0)
    }
    return(list(mnewdraw = mnewdraw, positionnew = positionnew, 
        addi = addvar, dropi = dropvar))
}
.rev.jump.int <-
function (molddraw = molddraw, K = K, mPlus = mPlus, maxk = Inf, 
    oldk = 0) 
{
    rev.idx = floor(runif(1, 0, 1) * 2)
    if ((rev.idx) | oldk == 0) {
        birth.death = .fls.samp.int(molddraw = molddraw, K = K, 
            mPlus = mPlus, maxk, oldk)
        mnewdraw = birth.death$mnewdraw
        positionnew = birth.death$positionnew
        addvar = birth.death$addi
        dropvar = birth.death$dropi
    }
    else {
        var.in = (1:K)[as.logical(molddraw)]
        var.out = (1:K)[!as.logical(molddraw)]
        mnewdraw = (molddraw > mPlus[, var.in[ceiling(length(var.in) * 
            runif(1, 0, 1))]])
        mnewdraw = mnewdraw | mPlus[var.out[ceiling(length(var.out) * 
            runif(1, 0, 1))], ]
        positionnew = (1:K)[mnewdraw]
        addvar = (1:K)[molddraw < mnewdraw]
        dropvar = (1:K)[molddraw > mnewdraw]
        if (length(dropvar) == 0) 
            dropvar = 0
        if (length(addvar) == 0) 
            addvar = 0
    }
    if (length(positionnew) > maxk) {
        return(.rev.jump.int(molddraw = molddraw, K = K, mPlus = mPlus, 
            maxk, oldk))
    }
    else {
        return(list(mnewdraw = as.numeric(mnewdraw), positionnew = positionnew, 
            addi = addvar, dropi = dropvar))
    }
}
.starter <-
function (K, start.value, y, N = N, XtX.big = XtX.big, Xty.big = Xty.big, 
    X = X, fixed.pos = numeric(0)) 
{
    if (is.na(start.value[1])) {
        start.value = min((N - 3), K)
    }
    if (any(start.value < 0) | !(is.numeric(start.value) | is.logical(start.value))) {
        start.value = min((N - 3), K)
        warning("Argument 'start.value' did not conform to required format. start.value has been changed to default - min(N-3,K)")
    }
    if (length(start.value) == 0) {
        start.value = numeric(K)
    }
    if (length(start.value) == 1) {
        if (start.value == 0) {
            start.value = numeric(K)
        }
    }
    if (length(start.value) > 1 && any(start.value > 1)) {
        sv = numeric(K)
        sv[start.value] = 1
        start.value = sv
        rm(sv)
    }
    if (length(start.value) == 1) {
        if (start.value > min((N - 3), K)) {
            cat("Submitted Start value is too large, used\n min(N-3,K) as starting model size instead\n\n")
            start.value = min((N - 3), K)
        }
        sorter = runif(K)
        start.position = order(sorter, seq(1:K))[1:start.value]
        XtX.start <- XtX.big[start.position, start.position]
        XtXinv.start <- chol2inv(chol(XtX.start))
        bhat = XtXinv.start %*% Xty.big[start.position]
        e = y - X[, start.position] %*% bhat
        sse = crossprod(e)
        s2 = as.numeric(sse/(N - length(start.position)))
        bcov = s2 * XtXinv.start
        bt = bhat/sqrt(diag(bcov))
        molddraw = rep(0, K)
        goodguy = as.numeric(abs(bt) > 0.2)
        molddraw[start.position] = goodguy
        start.position = (1:K)[as.logical(molddraw)]
        outstart = list(molddraw = molddraw, bhat = bhat, start.position = start.position)
    }
    if (length(start.value) > 1 && sum(start.value) == 0) {
        outstart = list(molddraw = rep(0, K), bhat = rep(0, K), 
            start.position = integer(0))
    }
    if (length(start.value) > 1 && sum(start.value) > 0) {
        if (length(start.value) != K) {
            stop("Starting Model contains unequal to K regressors,please respecify")
        }
        start.position = which(as.logical(start.value))
        XtX.start <- XtX.big[start.position, start.position]
        XtXinv.start <- chol2inv(chol(XtX.start))
        bhat = XtXinv.start %*% Xty.big[start.position]
        molddraw = rep(0, K)
        molddraw[start.position] = 1
        outstart = list(molddraw = molddraw, bhat = bhat, start.position = start.position)
    }
    fixed.pos = (1:K)[fixed.pos]
    if (length(fixed.pos) > 0) {
        outstart$molddraw[fixed.pos] = 1
        outstart$start.position = (1:K)[as.logical(outstart$molddraw)]
    }
    return(outstart)
}
.top10 <-
function (nmaxregressors = 10, nbmodels = 10, bbeta = FALSE, 
    lengthfixedvec = 0, bbeta2 = FALSE, ..., inivec_lik = numeric(0), 
    inivec_bool = character(0), inivec_count = numeric(0), inivec_vbeta = numeric(0), 
    inivec_vbeta2 = numeric(0), inivec_veck = 0, inivec_fixvec = numeric(0)) 
{
    findex = function() {
        seq_incl = seq_len(nbmodel)
        if (nbmodel == nbmodels) {
            seq_incl[indices] = seq_incl
        }
        else {
            truncindex = indices
            truncindex[(nbmodel + 1):nbmodels] = 0L
            seq_incl[truncindex] = seq_incl
        }
        return(seq_incl)
    }
    betamat = function(top10_betavec) {
        bins = (sapply(as.list(top10_bool[findex()]), hexobject$as.binvec))
        betamatx = matrix(0, nmaxregressors, nbmodel)
        if (length(top10_betavec) > 0) {
            betamatx[which(bins == 1)] = top10_betavec
        }
        else betamatx = betamatx[, 1]
        return(betamatx)
    }
    hexobject <- .hexcode.binvec.convert(nmaxregressors)
    if (nbmodels < 0) {
        nbmodels = 0
    }
    indices = integer(nbmodels)
    top10_lik = rep(-Inf, nbmodels)
    top10_bool = character(nbmodels)
    top10_count = integer(nbmodels)
    top10_fixvec = numeric(lengthfixedvec * nbmodels)
    if (bbeta) 
        lbetas = vector("list", nbmodels)
    if (bbeta2) 
        lbetas2 = vector("list", nbmodels)
    seq_nbmodel = seq_len(nbmodels)
    ix_of_mybool = logical(nbmodels)
    nbmodel = length(inivec_lik)
    top10_lik[seq_len(nbmodel)] = inivec_lik
    top10_count[seq_len(nbmodel)] = inivec_count
    if (is.character(inivec_bool)) {
        top10_bool[seq_len(nbmodel)] = inivec_bool
    }
    else {
        if (is.vector(inivec_bool) & (length(inivec_bool) == 
            nmaxregressors)) {
            top10_bool[seq_len(nbmodel)] = hexobject$as.hexcode(inivec_bool)
        }
        else if (is.list(inivec_bool)) {
            top10_bool[seq_len(nbmodel)] = sapply(inivec_bool, 
                hexobject$as.hexcode)
        }
        else if (is.matrix(inivec_bool)) {
            top10_bool[seq_len(nbmodel)] = sapply(as.list(as.data.frame(inivec_bool)), 
                hexobject$as.hexcode)
        }
        else stop("inivec_bool is wrong format!")
    }
    top10_fixvec = inivec_fixvec
    if (is.na(inivec_veck[1])) {
        inivec_veck = 0
    }
    if (bbeta | bbeta2) {
        veck_ix = c(0, cumsum(inivec_veck))
        veckix_aux = as.list(seq_len(nbmodel))
        veckix_aux = lapply(veckix_aux, function(x) {
            if (veck_ix[[x]] == veck_ix[[x + 1]]) 
                c(0, 0)
            else c(veck_ix[[x]] + 1, veck_ix[[x + 1]])
        })
    }
    if (bbeta) {
        lbetas[seq_len(nbmodel)] = lapply(veckix_aux, function(x) inivec_vbeta[x[[1]]:x[[2]]])
    }
    else lbetas = list(numeric(0))
    if (bbeta2) {
        lbetas2[seq_len(nbmodel)] = lapply(veckix_aux, function(x) inivec_vbeta2[x[[1]]:x[[2]]])
    }
    else lbetas2 = list(numeric(0))
    lastvec01 = integer(nmaxregressors)
    modidx = length(top10_lik)
    indices[seq_len(nbmodel)] = order(inivec_lik, decreasing = TRUE)
    min.index = which.max(indices)
    if (length(min.index) > 0) {
        min.top10_lik = top10_lik[[min.index]]
    }
    else {
        if (nbmodels > 0) 
            min.top10_lik = -Inf
        else min.top10_lik = Inf
    }
    index.of.mybool = function(mybool) {
        ix_of_mybool <<- (mybool == top10_bool)
    }
    check4dupl = index.of.mybool
    dupl.possible = TRUE
    retlist = list(addmodel = function(mylik, vec01, vbeta = numeric(0), 
        vbeta2 = numeric(0), fixedvec = numeric(0)) {
        if (mylik >= min.top10_lik | nbmodel < nbmodels) {
            if (identical(lastvec01, vec01)) {
                top10_count[[modidx]] <<- top10_count[[modidx]] + 
                  1
            } else {
                lastvec01 <<- vec01
                mybool = hexobject$as.hexcode(vec01)
                check4dupl(mybool)
                if (!any(ix_of_mybool)) {
                  if (nbmodel < nbmodels) {
                    nbmodel <<- nbmodel + 1
                    modidx <<- nbmodel
                  } else {
                    modidx <<- min.index
                  }
                  ltmylik = (top10_lik <= mylik)
                  indices <<- indices + ltmylik
                  indices[[modidx]] <<- nbmodels - sum(ltmylik) + 
                    1
                  top10_lik[[modidx]] <<- mylik
                  top10_bool[[modidx]] <<- mybool
                  top10_count[[modidx]] <<- 1
                  min.index <<- which.max(indices)
                  min.top10_lik <<- top10_lik[[min.index]]
                  if (lengthfixedvec > 0) {
                    top10_fixvec[(modidx - 1) * lengthfixedvec + 
                      seq_len(lengthfixedvec)] <<- fixedvec
                  }
                  if (bbeta) {
                    lbetas[[modidx]] <<- vbeta
                  }
                  if (bbeta2) {
                    lbetas2[[modidx]] <<- vbeta2
                  }
                } else {
                  modidx <<- seq_nbmodel[ix_of_mybool]
                  top10_count[[modidx]] <<- top10_count[[modidx]] + 
                    1
                }
            }
        }
    }, lik = function() {
        return(top10_lik[findex()])
    }, bool = function() {
        return(top10_bool[findex()])
    }, ncount = function() {
        return(top10_count[findex()])
    }, nbmodels = nbmodels, nregs = nmaxregressors, betas_raw = function() {
        return(unlist(lbetas[findex()]))
    }, betas2_raw = function() {
        return(unlist(lbetas2[findex()]))
    }, kvec_raw = function() {
        return(sapply(lbetas, length)[findex()])
    }, bool_binary = function() {
        return(sapply(as.list(top10_bool[findex()]), hexobject$as.binvec))
    }, betas = function() {
        betamat(unlist(lbetas[findex()]))
    }, betas2 = function() {
        betamat(unlist(lbetas2[findex()]))
    }, fixed_vector = function() {
        if (lengthfixedvec <= 0) {
            return(matrix(0, 0, 0))
        } else {
            findex_base = (findex() - 1) * lengthfixedvec
            findex_fixvec = numeric(0)
            for (xx in 1:lengthfixedvec) findex_fixvec = rbind(findex_fixvec, 
                findex_base + xx)
            return(matrix(top10_fixvec[c(findex_fixvec)], lengthfixedvec))
        }
    }, duplicates_possible = function(possible = NULL) {
        if (!is.logical(possible)) return(dupl.possible)
        if (possible) {
            check4dupl <<- index.of.mybool
            dupl.possible <<- TRUE
            ix_of_mybool <<- logical(nbmodels)
        } else {
            check4dupl <<- function(mybool) {
            }
            dupl.possible <<- FALSE
            ix_of_mybool <<- FALSE
        }
    })
    class(retlist) = "topmod"
    return(retlist)
}
.topmod.as.bbetaT <-
function (tm, gprior.info = NULL, yXdata = NULL, addr2 = FALSE) 
{
    is.bmao = FALSE
    if (is.bma(tm)) {
        is.bmao = TRUE
        bmao = tm
        yXdata = bmao$X.data
        gprior.info = bmao$gprior.info
        tm = bmao$topmod
    }
    yXdata = as.matrix(yXdata)
    N = nrow(yXdata)
    K = ncol(yXdata) - 1
    yXdata = yXdata - matrix(colMeans(yXdata), N, K + 1, byrow = TRUE)
    if (length(tm$lik()) < 1) {
        if (is.bmao) 
            return(bmao)
        else return(tm)
    }
    if (!addr2) 
        if ((length(tm$betas_raw()) > 0) & (ncol(as.matrix(tm$betas())) == 
            length(tm$lik()))) {
            if (is.bmao) 
                return(bmao)
            else return(tm)
        }
    bools = (tm$bool_binary())
    yty = c(crossprod(yXdata[, 1]))
    positions = lapply(lapply(as.list(as.data.frame(bools)), 
        as.logical), which)
    olsmodels = lapply(lapply(positions, .ols.terms2, yty = yty, 
        N = N, K = K, XtX.big = crossprod(yXdata[, -1]), Xty.big = c(crossprod(yXdata[, 
            -1], yXdata[, 1]))), function(x) x$full.results())
    lprobo = gprior.info$lprobcalc
    lpl = lapply(olsmodels, function(x) lprobo$lprob(x$ymy, length(x$bhat), 
        x$bhat, x$diag.inverse))
    veck = as.vector(unlist(lapply(lapply(lpl, "[[", "b1new"), 
        length)))
    b1raw = as.vector(unlist(lapply(lpl, "[[", "b1new")))
    b2raw = as.vector(unlist(lapply(lpl, "[[", "b2new")))
    fixedvecmat = tm$fixed_vector()
    if (addr2) {
        r2 = 1 - sapply(olsmodels, function(x) x$ymy)/yty
        if (nrow(fixedvecmat) == 0) {
            fixedvecmat = matrix(0, 0, length(veck))
        }
        else if (mean(abs(r2 - fixedvecmat[1, ])) < 1e-17) {
            fixedvecmat = fixedvecmat[-1, , drop = FALSE]
        }
        fixedvecmat = rbind(r2, fixedvecmat)
    }
    lengthfixedvec = nrow(fixedvecmat)
    tm <- .top10(nmaxregressors = tm$nregs, nbmodels = tm$nbmodels, 
        bbeta = TRUE, lengthfixedvec = lengthfixedvec, bbeta2 = TRUE, 
        inivec_lik = tm$lik(), inivec_bool = tm$bool(), inivec_count = tm$ncount(), 
        inivec_vbeta = b1raw, inivec_vbeta2 = b2raw, inivec_veck = veck, 
        inivec_fixvec = c(fixedvecmat))
    if (is.bmao) {
        bmao$topmod <- tm
        return(bmao)
    }
    return(tm)
}
