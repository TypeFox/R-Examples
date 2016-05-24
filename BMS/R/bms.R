bms <-
function (X.data = NULL, burn = 1000, iter = NA, nmodel = 500, 
    mcmc = "bd", g = "UIP", mprior = "random", mprior.size = NA, 
    user.int = TRUE, start.value = NA, g.stats = TRUE, logfile = FALSE, 
    logstep = 10000, force.full.ols = FALSE, fixed.reg = numeric(0), 
    data = NULL, randomizeTimer = TRUE) 
{
    if (missing(X.data) & !missing(data)) 
        X.data = data
    mf <- match.call(expand.dots = FALSE)
    if (!is.na(match("X.data", names(mf)))) 
        names(mf)[[match("X.data", names(mf))]] = "formula"
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    if (is.data.frame(X.data)) {
        mf <- X.data
    }
    else if (is.matrix(X.data)) {
        mf <- model.frame(as.data.frame(X.data, drop.unused.levels = TRUE))
    }
    else {
        mf <- eval(mf, parent.frame())
    }
    X.data = as.matrix(mf)
    if (any(is.na(X.data))) {
        X.data = na.omit(X.data)
        if (nrow(X.data) < 3) {
            stop("Too few data observations. Please provide at least three data rows without NA entries.")
        }
        warning("Argument 'X.data' contains NAs. The corresponding rows have not been taken into account.")
    }
    N <- nrow(X.data)
    K = ncol(X.data) - 1
    maxk = N - 3
    if (is.null(nmodel[1]) || is.na(nmodel[1]) || nmodel[1] <= 
        0) {
        dotop = FALSE
        nmodel = 0
    }
    else {
        dotop = TRUE
    }
    nameix = 1:K
    names(nameix) = colnames(X.data[, -1, drop = FALSE])
    fixed.pos = nameix[fixed.reg]
    rm(nameix)
    if (missing(mcmc) && ((K - length(fixed.pos)) < 15)) {
        mcmc = "enum"
    }
    int = FALSE
    is.enum = FALSE
    if (is.function(mcmc)) {
        samplingfun = mcmc
        mcmc = "custom"
    }
    else {
        if (length(grep("int", mcmc, ignore.case = TRUE))) {
            int = TRUE
        }
        if (length(grep("enum", mcmc, ignore.case = TRUE))) {
            is.enum = TRUE
            samplingfun = .iterenum
            if (K > maxk) 
                samplingfun = .iterenum.KgtN
        }
        else if (length(grep("bd", mcmc, ignore.case = TRUE))) {
            samplingfun = switch(int + 1, .fls.samp, .fls.samp.int)
        }
        else {
            samplingfun = switch(int + 1, .rev.jump, .rev.jump.int)
        }
    }
    if (int && (length(fixed.pos) > 0L)) {
        warning("interaction sampler does not allow for non-zero argument fixed.pos; consequently it was set fixed.pos=0")
        fixed.pos = numeric(0)
    }
    sampling = .fixedset.sampler(samplingfun, fullK = K, fixed.pos = fixed.pos, 
        X.data = X.data)
    if (is.enum) {
        burn = 0
        int = FALSE
        mcmc = "enum"
        is.enum = TRUE
        tmp = .enum_startend(iter = iter, start.value = start.value, 
            K = K, maxk = maxk, fixed.pos = fixed.pos)
        iter = tmp$iter
        start.value = tmp$start.value
    }
    else {
        if (is.na(iter)) {
            iter = 3000
        }
    }
    if (logfile != FALSE) {
        if (is.character(logfile)) {
            sfilename = logfile
        }
        else {
            sfilename = ""
        }
        if (nchar(sfilename) > 0) 
            file.create(sfilename)
        logfile = TRUE
        cat(as.character(Sys.time()), ": starting loop ... \n", 
            append = TRUE, file = sfilename)
        if (logstep != 10000) 
            fact = logstep
        else fact = max(floor((burn + iter)/100), logstep)
    }
    y <- as.matrix(X.data[, 1])
    X <- as.matrix(X.data[, 2:ncol(X.data)])
    y <- y - mean(y)
    X <- X - matrix(colMeans(X), N, K, byrow = TRUE)
    XtX.big = crossprod(X)
    Xty.big = crossprod(X, y)
    yty = as.vector(crossprod(y))
    coreig = eigen(cor(X), symmetric = TRUE, only.values = TRUE)$values
    if (!force.full.ols) {
        if (sum(coreig > 1e-07) < min(K, (N - 1))) {
            force.full.ols = TRUE
        }
    }
    if (any(coreig[1:min(K, (N - 1))] < 1e-16)) {
        warning(paste("data seems to be rank-deficient: its rank seems to be only ", 
            sum(coreig > 1e-13)))
    }
    if (int) {
        if (length(grep("#", colnames(X.data), fixed = TRUE)) == 
            0) 
            stop("Please separate column names of interaction terms by # (e.g. A#B)")
        mPlus = .constr.intmat(X, K)
    }
    else {
        mPlus <- NA
    }
    pmplist = .choose.mprior(mprior, mprior.size, K = K, X.data = X.data, 
        fixed.pos = fixed.pos)
    mprior = pmplist$mp.mode
    gprior.info = .choose.gprior(g, N, K, return.g.stats = g.stats, 
        yty = yty, X.data = X.data)
    lprobcalc = gprior.info$lprobcalc
    start.list = .starter(K, start.value, y, N = N, XtX.big = XtX.big, 
        Xty.big = Xty.big, X = X, fixed.pos = fixed.pos)
    molddraw = start.list$molddraw
    start.position = start.list$start.position
    kold = sum(molddraw)
    position = (1:K)[molddraw == 1]
    collect.otherstats = FALSE
    otherstats = numeric(0)
    add.otherstats = numeric(0)
    if (gprior.info$return.g.stats & !(gprior.info$is.constant)) {
        add.otherstats = gprior.info$shrinkage.moments
        collect.otherstats = TRUE
    }
    cumsumweights = iter
    null.lik = lprobcalc$just.loglik(ymy = yty, k = 0)
    if (collect.otherstats) {
        addup <- function() {
            inccount <<- inccount + molddraw
            msize <<- msize + kold
            if (kold != 0) {
                bm[c(position, K + position, 2 * K + kold, 3 * 
                  K + position)] = c(b1, b2, 1, b1 > 0)
                bmo <<- bmo + bm
            }
            else {
                null.count <<- null.count + 1
            }
            otherstats <<- lik.list[["otherstats"]]
            add.otherstats <<- add.otherstats + otherstats
        }
    }
    else {
        addup <- function() {
            inccount <<- inccount + molddraw
            msize <<- msize + kold
            if (kold != 0) {
                bm[c(position, K + position, 2 * K + kold, 3 * 
                  K + position)] = c(b1, b2, 1, b1 > 0)
                bmo <<- bmo + bm
            }
            else {
                null.count <<- null.count + 1
            }
        }
    }
    if (is.enum) {
        cumsumweights = 0
        if (collect.otherstats) {
            addup <- function() {
                weight = exp(pmpold + lprobold - null.lik)
                inccount <<- inccount + weight * molddraw
                msize <<- msize + weight * kold
                cumsumweights <<- cumsumweights + weight
                if (kold != 0) {
                  bm[c(position, K + position, 2 * K + kold, 
                    3 * K + position)] = weight * c(b1, b2, 1, 
                    b1 > 0)
                  bmo <<- bmo + bm
                }
                else {
                  null.count <<- null.count + weight
                }
                otherstats <<- lik.list[["otherstats"]]
                add.otherstats <<- add.otherstats + weight * 
                  otherstats
            }
        }
        else {
            addup <- function() {
                weight = exp(pmpold + lprobold - null.lik)
                inccount <<- inccount + weight * molddraw
                msize <<- msize + weight * kold
                cumsumweights <<- cumsumweights + weight
                if (kold != 0) {
                  bm[c(position, K + position, 2 * K + kold, 
                    3 * K + position)] = weight * c(b1, b2, 1, 
                    b1 > 0)
                  bmo <<- bmo + bm
                }
                else {
                  null.count <<- null.count + weight
                }
            }
        }
    }
    environment(addup) <- environment()
    ols.object = .ols.terms2(positions = (1:K)[molddraw == 1], 
        yty = yty, k = kold, N, K = K, XtX.big = XtX.big, Xty.big = Xty.big)
    lik.list = lprobcalc$lprob.all(ymy = ols.object$ymy, k = kold, 
        bhat = ols.object$bhat, diag.inverse = ols.object$diag.inverse)
    lprobold = lik.list$lprob
    b1 = lik.list$b1new
    b2 = lik.list$b2new
    pmpold = pmplist$pmp(ki = kold, mdraw = molddraw)
    topmods = topmod(nbmodels = nmodel, nmaxregressors = K, bbeta = FALSE, 
        lengthfixedvec = length(add.otherstats))
    if (mcmc == "enum") {
        try(topmods$duplicates_possible(FALSE), silent = TRUE)
    }
    if (dotop && (burn == 0L)) 
        topmods$addmodel(mylik = pmpold + lprobold, vec01 = molddraw, 
            fixedvec = lik.list$otherstats)
    null.count = 0
    models.visited = 0
    inccount = numeric(K)
    msize = 0
    k.vec = numeric(K)
    b1mo = numeric(K)
    ab = numeric(K)
    b2mo = numeric(K)
    bb = numeric(K)
    possign = inccount
    mnewdraw = numeric(K)
    if (force.full.ols) {
        candi.is.full.object = TRUE
    }
    else {
        candi.is.full.object = FALSE
    }
    bmo = numeric(4 * K)
    bm = bmo
    if (is.enum) {
        addup()
    }
    if (!is.finite(pmpold)) 
        pmpold = -1e+90
    if (randomizeTimer) 
        set.seed(as.numeric(Sys.time()))
    t1 <- Sys.time()
    nrep = burn + iter
    i = 0
    while (i < nrep) {
        i = i + 1
        if (logfile) {
            if (i%%fact == 0) {
                if (nchar(sfilename) == 0) 
                  message(as.character(Sys.time()), ":", i, "current draw")
                else cat(as.character(Sys.time()), ":", i, "current draw \n", 
                  append = TRUE, file = sfilename)
            }
        }
        a = sampling(molddraw = molddraw, K = K, mPlus = mPlus, 
            maxk = maxk, oldk = kold)
        mnewdraw = a[["mnewdraw"]]
        positionnew = a[["positionnew"]]
        knew = length(positionnew)
        pmpnew = pmplist[["pmp"]](ki = knew, mdraw = mnewdraw)
        if (!is.enum) {
            if (int) {
                if (length(c(a$dropi, a$addi)) > 2 | i < 3 | 
                  force.full.ols) {
                  candi.is.full.object = TRUE
                }
                else {
                  candi.is.full.object = FALSE
                }
            }
            if (candi.is.full.object) {
                ols.candidate = .ols.terms2(positions = positionnew, 
                  yty = yty, k = knew, N, K = K, XtX.big = XtX.big, 
                  Xty.big = Xty.big)
                ymy.candi = ols.candidate[["ymy"]]
            }
            else {
                ymy.candi = ols.object[["child.ymy"]](a$addi, 
                  a$dropi, k = knew)
            }
            if ((ymy.candi < 0) | is.na(ymy.candi)) 
                stop(paste("stumbled on rank-deficient model"))
            lprobnew = lprobcalc[["just.loglik"]](ymy = ymy.candi, 
                k = knew)
            accept.candi = as.logical(log(runif(1, 0, 1)) < lprobnew - 
                lprobold + pmpnew - pmpold)
        }
        else {
            accept.candi = TRUE
            candi.is.full.object = FALSE
        }
        if (accept.candi) {
            if (!candi.is.full.object) {
                ols.res = ols.object[["mutate"]](addix = a$addi, 
                  dropix = a$dropi, newpos = positionnew, newk = knew)
            }
            else {
                ols.object = ols.candidate
                ols.res = ols.candidate[["full.results"]]()
            }
            lik.list = lprobcalc[["lprob.all"]](max(0, ols.res$ymy), 
                knew, ols.res$bhat, ols.res$diag.inverse)
            lprobold = lik.list[["lprob"]]
            position = positionnew
            pmpold = pmpnew
            molddraw = mnewdraw
            kold = knew
            models.visited = models.visited + 1
        }
        if (i > burn) {
            b1 = lik.list[["b1new"]]
            b2 = lik.list[["b2new"]]
            addup()
            if (dotop) 
                topmods[["addmodel"]](mylik = pmpold + lprobold, 
                  vec01 = molddraw, fixedvec = otherstats)
        }
    }
    if (dotop) 
        topmods = .topmod.as.bbetaT(topmods, gprior.info, X.data)
    timed <- difftime(Sys.time(), t1)
    if (is.enum) {
        iter = iter + 1
        models.visited = models.visited + 1
    }
    bmo = matrix(bmo, 4, byrow = TRUE)
    b1mo = bmo[1, ]
    b2mo = bmo[2, ]
    k.vec = bmo[3, ]
    possign = bmo[4, ]
    rm(bmo)
    post.inf = .post.calc(gprior.info, add.otherstats, k.vec, 
        null.count, X.data, topmods, b1mo, b2mo, iter, burn, 
        inccount, models.visited, K, N, msize, timed, cumsumweights, 
        mcmc, possign)
    result = list(info = post.inf$info, arguments = .construct.arglist(bms, 
        environment()), topmod = topmods, start.pos = sort(start.position), 
        gprior.info = post.inf$gprior.info, mprior.info = pmplist, 
        X.data = X.data, reg.names = post.inf$reg.names, bms.call = try(match.call(bms, 
            sys.call(0)), silent = TRUE))
    class(result) = c("bma")
    if (user.int) {
        print(result)
        print(timed)
        plot.bma(result)
    }
    return(invisible(result))
}
