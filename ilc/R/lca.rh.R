lca.rh <-
function(dat, year=dat$year, age=dat$age, series=1, max.age=100, 
                   dec.conv=6, clip = 3, error=c('poisson', 'gaussian'), 
                   model = c('m', 'h0', 'h1', 'h2', 'ac', 'lc'),  
                   restype = c("logrates", "rates", "deaths", "deviance"), 
                   scale = F, interpolate = F, verbose=T, spar = NULL){
    
    # --- Preamble: --- #
    if (class(dat) != "demogdata" || dat$type != "mortality")
        stop("Not mortality data")
    error.choices <- c('poisson', 'gaussian')
    if (is.character(error)) error <- match.arg(error, error.choices)
    else error <- error.choices[error]
    res.choices <- c("logrates", "rates", "deaths", "deviance")
    if (is.character(restype)) {
        restype <- match.arg(restype, res.choices)
        restype <- match(restype, res.choices)
    } else {
        restype <- as.integer(restype)
        if (restype < 1 || restype > 4) stop('\n  Unspecified residual type! ')
    }
    mod.choices <- c('m', 'h0', 'h1', 'h2', 'ac', 'lc')
    if (is.character(model)) { 
        model <- match.arg(model, mod.choices)
        model <- match(model, mod.choices)
    } else {
        model <- as.integer(model)
        if (model < 1 || model > 6) stop('\n  Unspecified model! ')
    }
    
    # setup:
    dum = lca.set(dat, year, age, series, max.age, interpolate)
    
    n = dum$n
    k = dum$k
    rtx = dum$rtx
    mx = dum$mx
    md = dum$md
    mz = dum$mz
    me = dum$me
    
    ft <- gl(n, k) 
    if (model <= 5){
        stx <- year[1]-age[k]   # year of birth of first cohort  
        etx <- year[n]-age[1]   # year of birth of last cohort 
        mtx <- array(year[ft]-age, dim(mx)) 
        fx <- factor(rep(seq(age), n))
        ftx <- factor(mtx - stx + 1)
    }
    
    mod.choices <- c(' M = a(x)+b0(x)*i(t-x)+b1(x)*k(t) ', ' H(0) = a(x)+i(t-x)+k(t) ',
                     ' H(1) = a(x)+i(t-x)+b1(x)*k(t) ', ' H(2) = a(x)+b0(x)*i(t-x)+k(t) ', 
                     ' AC = a(x)+b0(x)*i(t-x) ', ' LC = a(x)+b1(x)*k(t) ')    
    cat('\n  Fitting model:', sqb(mod.choices[model]), '\n\t- with', error, 
        'error structure', if (error=='poisson') 'and with deaths as weights', 
        '-\n')
    
    if (verbose){
        cat('Note:', sum(bool(md < 1e-9, na=T)), 'cells have 0/NA deaths and ',
            sum(!mw), 'have 0/NA exposure', '\n  out of a total of', n*k, 
            'data cells.\n')
    }
    
    # Clip the 'corner' cohorts:
    if (clip && model <= 5){ 
        mw[bool(mtx < stx+clip)] <- F   
        mw[bool(mtx > etx-clip)] <- F
        warning(paste('The cohorts outside',  sqb(paste(c(stx+clip, etx-clip), 
                                                        collapse=', ')), 'were zero weighted (clipped).'))
    }
    
    # --- End Preamble. Start fitting: --- #
    
    # Get starting values:
    # this gives the log of the geometrical means of the mortality rates:
    # note that this is a weighted mean (with the exception of the LC model):
    ax <- mz;  if (model <= 5) ax <- ax*mw
    ax <- apply(ax, 1, mean, na.rm=T) 
    
    # fitting the model as a single vector using the defined factors (see setup)
    kt <- rep(0, n); names(kt) <- year
    itx <- rep(0, rtx); names(itx) <- seq(year[1]-age[k], l=rtx)
    if (model < 5) {
        if (verbose) cat(' Automated start: initial values by glm fitting',  
                         'of factors i(t-x)+k(t).')
        if (error == 'poisson'){
            y <- md; dim(y) <- NULL
            os <- log(me)+ax; dim(os) <- NULL
            w <- as.numeric(mw)
            fit <- glm(y ~ offset(os)+ftx+ft-1, family=poisson, weights=w) # family=error
        }  else {
            y <- mz; dim(y) <- NULL
            fit <- glm(y ~ ftx+ft-1, family=gaussian)  # fam=error
        }
        fit.coef <- names(fit$coef)
        ind <- grep('ftx[0-9]', fit.coef)
        itx[as.numeric(substring(fit.coef[ind], 4))] <- fit$coef[ind]
        itx[is.na(itx)] <- 0
        ind <- grep('ft[0-9]', fit.coef)
        kt[as.numeric(substring(fit.coef[ind], 3))] <- fit$coef[ind]
        kt[is.na(kt)] <- 0
    }
    if (model == 5) {
        if (verbose) cat(' Automated start: initial values by glm fitting',  
                         'of factors ax+i(t-x).')
        if (error == 'poisson'){
            y <- md; dim(y) <- NULL
            os <- log(me); dim(os) <- NULL
            w <- as.numeric(mw)
            fit <- glm(y ~ offset(os)+fx+ftx-1, family=poisson, weights=w) # fam=error
        }  else {
            y <- mz; dim(y) <- NULL
            w <- mw*bool(md > 1e-9, na=F); dim(w) <- NULL
            fit <- glm(y ~ fx+ftx-1, family=gaussian, weights=w)  # fam=error
        }
        fit.coef <- names(fit$coef)
        ind <- grep('fx[0-9]', fit.coef)
        ax[as.numeric(substring(fit.coef[ind], 3))] <- fit$coef[ind]
        ax[is.na(ax)] <- 0
        ind <- grep('ftx[0-9]', fit.coef)
        itx[as.numeric(substring(fit.coef[ind], 4))] <- fit$coef[ind]
        itx[is.na(itx)] <- 0
    }
    
    # set bx0 and bx1 starting values depending on models:
    bx0 <- switch(model, rep(1, k), rep(1, k), rep(1, k),       # m h0 h1 
                  rep(1/k, k), rep(1/k, k), rep(0, k))  # h2 ac lc
    bx1 <- switch(model, rep(1, k), rep(1, k), rep(1/k, k),     # m h0 h1 
                  rep(1, k), rep(0, k), rep(1/k, k))    # h2 ac lc  
    names(bx0) <- names(bx1) <- age
    if (verbose) {
        cat('\n Starting values are:\n')
        mcoef <- list(age=ax, bx1=bx1, per=kt)
        if (model <= 5) {
            i <- round(rtx/max(n,k))
            if (i >= 2){
                cc <- split(seq(rtx), cut(seq(rtx), i, order=T))
                names(cc) <- paste('coh', seq(cc), sep='')
                for(j in seq(i)) cc[[j]] <- itx[cc[[j]]]
            } else cc <- list(coh=itx)
            mcoef$bx0 <- bx0
            mcoef <- c(mcoef, cc)
        }
        class(mcoef) <- 'coef'
        print(mcoef)
    }
    
    # Assign target/observed values (logrates or number of deaths):
    if (error == 'gaussian') { my <- mz } 
    if (error == 'poisson') { my <- md } 
    # fit LCM (Gaussian error structure, central rates as target)
    mf <- lcc(ax, bx0, itx, bx1, kt)
    # fit LCM (Poisson error structure, deaths as target)
    if (error == 'poisson') mf <- me*exp(mf)
    # Calculate deviance:
    D <- deviance.lca(my, mf, mw, error)    
    if (verbose) cat('\n Iterative fit:\n #iter   Dev    non-conv\n')
    dD <- 1 # a fictive difference of 1
    ncv <- iter <- 0  # non-convergence and iteration counters
    while(round(dD, dec.conv)){   
        # Note: with higher precisions it might over-fit the parameters 
        iter <- iter + 1
        if (verbose) cat('  ', c(iter, round(D, 6), ncv), '\n ', sep='  ')
        # --------  Stage 0: update ax ----------
        if (model >= 5){
            td <- apply(mw*(my-mf), 1, sum, na.rm=T)
            edt <- mw
            if (error == 'poisson') edt <- edt*mf
            edt <- apply(edt, 1, sum, na.rm=T)
            ax <- ax + td/edt
            # fit LC (Gaussian error structure, central rates as target)
            mf <- lcc(ax, bx0, itx, bx1, kt)
            # fit LC (Poisson error structure, deaths as target)
            if (error == 'poisson') mf <- me*exp(mf)
        }
        # --------  Stage 1: update itx ----------
        if (model <= 5){
            td <- mw*(my-mf)*bx0
            edt <- mw*bx0^2
            if (error == 'poisson') edt <- edt*mf
            for(i in seq(rtx)){
                ind <- ftx==i
                u.itx <- sum(edt[ind], na.rm=T)
                if (u.itx) u.itx <- sum(td[ind], na.rm=T)/u.itx
                itx[i] <- itx[i] + u.itx
            }
            if (model == 5) itx <- itx - itx[1]  
            # fit LCM (Gaussian error structure, central rates as target)
            mf <- lcc(ax, bx0, itx, bx1, kt)
            # fit LCM (Poisson error structure, deaths as target)
            if (error == 'poisson') mf <- me*exp(mf)
        }
        # --------  Stage 2: update bx0 ----------
        if (model == 1 || model == 4 || model == 5){
            td <- apply(mw*(my-mf)*itx[ftx], 1, sum, na.rm=T)
            edt <- mw*itx[ftx]^2
            if (error == 'poisson') edt <- edt*mf
            edt <- apply(edt, 1, sum, na.rm=T)
            bx0 <- bx0 + td/edt
            if (!is.null(spar)) bx0 <- smooth.spline(bx0, spar=spar)$y
            # fit LCM (Gaussian error structure, central rates as target)
            mf <- lcc(ax, bx0, itx, bx1, kt)
            # fit LCM (Poisson error structure, deaths as target)
            if (error == 'poisson') mf <- me*exp(mf)
        }
        # --------  Stage 3: update kt ----------
        if (model != 5){
            td <- apply(mw*(my-mf)*bx1, 2, sum, na.rm=T)
            edt <- mw*bx1^2
            if (error == 'poisson') edt <- edt*mf
            edt <- apply(edt, 2, sum, na.rm=T)
            kt <- kt + td/edt; 
            # adjust updates:
            adj <- if (model == 6) mean(kt) else kt[1]
            kt <-  kt-adj 
            # fit LCM (Gaussian error structure, central rates as target)
            mf <- lcc(ax, bx0, itx, bx1, kt)
            # fit LCM (Poisson error structure, deaths as target)
            if (error == 'poisson') mf <- me*exp(mf)
        }
        # --------  Stage 4: update bx1 ----------
        if (model == 1 || model == 3 || model == 6){
            td <- apply(mw*(my-mf)*kt[ft], 1, sum, na.rm=T)
            edt <- mw*kt[ft]^2
            if (error == 'poisson') edt <- edt*mf
            edt <- apply(edt,  1, sum, na.rm=T)
            bx1 <- bx1 + td/edt
            if (!is.null(spar)) bx1 <- smooth.spline(bx1, spar=spar)$y
            # fit LCM (Gaussian error structure, central rates as target)
            mf <- lcc(ax, bx0, itx, bx1, kt)
            # fit LCM (Poisson error structure, deaths as target)
            if (error == 'poisson') mf <- me*exp(mf)
        }
        # --------  Stage 5: check deviance convergence ----------
        dD <- D - deviance.lca(my, mf, mw, error)
        # Make the iteration non-convergent if the deviance increases 5 times 
        # in a row with more than 10 units (only a subjective criteria).
        if (dD < -10){
            if (ncv > 10) stop('Non-convergent deviance!\n')
            else ncv <- ncv+1
        } else ncv <- 0
        D <- D - dD
    }
    
    # reset the names of bx0 and bx1 if smoothing applied:
    if (!is.null(spar)) names(bx0) <- names(bx1) <- age
    
    cat('\n Iterations finished in:', iter, 'steps\n')
    
    if (model == 1 || model == 4 || model == 5){
        sb <- sum(bx0)
        bx0 <- bx0/sb; itx <- itx*sb
    }
    if (model == 1 || model == 3 || model == 6){
        sb <- sum(bx1)
        bx1 <- bx1/sb; kt <- kt*sb
    }
    
    # Prepare return object:
    mf <- lcc(ax, bx0, itx, bx1, kt)  # fitted logrates
    mdf <- me*exp(mf)     # fitted number of deaths
    # for deviance residuals:
    mD <- deviance.lca(my, if (error=='poisson') mdf else mf, mw, error, total=F)
    # degrees of freedom depending on model and adjusted for missing data cells:
    nu <- switch(model, k*(n-3)-(2*n-3), k*(n-1)-(2*n-1), k*(n-2)-2*(n-1),
                 k*(n-2)-2*(n-1), (k-1)*(n-3), (k-1)*(n-2)) - sum(!mw)   
    sc.phi <- sum(mD, na.rm=T)/nu  # scaling factor (phi)
    fit <- switch(restype, mf, exp(mf), mdf, if (error=='poisson') mdf else mf)
    fitted <- fts(age, fit, frequency = 1, start = year[1], xname = "Age", 
                  yname = paste("Fitted", 
                                rb(res.choices[if (restype<4) restype else if (error=='poisson') 3 else 1]), 
                                "mortality rate"))
    fitted$type <- res.choices[if (restype<4) restype else if (error=='poisson') 3 else 1]
    res <- switch(restype, mz, mx, md, my) - fit
    if (restype == 4) res <- sign(res)*sqrt(mD/sc.phi)    
    residuals <- fts(age, res, frequency = 1, start = year[1], xname = "Age",  
                     yname = paste("Residuals", rb(res.choices[restype]), "mortality rate"))
    residuals$type <- res.choices[restype]
    if (scale) {
        avdiffk <- -mean(diff(kt))
        bx1 <- bx1 * avdiffk
        kt <- kt/avdiffk
        # re-applied the above to cohorts factor 
        avdiffi <- -mean(diff(itx))
        bx0 <- bx0 * avdiffi
        itx <- itx/avdiffi
        warning('\"kt\", \"itx\", \"bx0\" and \"bx1\" parameters were re-scaled!')
    }
    if (verbose) {
        cat('\n Updated values are:\n')
        mcoef <- list(age=ax, bx1=bx1, per=kt)
        if (model <= 5) {
            i <- round(rtx/max(n,k))
            if (i >= 2){        
                cc <- split(seq(rtx), cut(seq(rtx), i, order=T))
                names(cc) <- paste('coh', seq(cc), sep='')
                for(j in seq(i)) cc[[j]] <- itx[cc[[j]]]
            } else cc <- list(coh=itx)
            mcoef$bx0 <- bx0
            mcoef <- c(mcoef, cc)
        }
        class(mcoef) <- 'coef'
        print(mcoef, dec=5)
        # print out total sums too:
        cat('\t total sums are: \n')                
        psu <- list(b0=bx0, b1=bx1, itx=itx, kt=kt)
        psu <- unlist(lapply(psu, sum))
        print(round(psu, 4))
    }
    
    # Alternative deviance scaling by Hyndman (see lca in demography package): 
    drift <- mean(diff(kt))
    kt.lin <- mean(kt) + drift * (1:n - (n + 1)/2)
    mdf.lin <- me * exp(lcc(ax, bx0, itx, bx1, kt.lin))
    nu.lin <- k * (n - 2)   # k/(k-1)*nu
    sc.phi.lin <- deviance.lca(md, mdf.lin, mw,'poisson')/nu.lin
    output <- list(label = dat$label, age = age, year = year, mx = mx, ax = ax, 
                   bx = bx1, kt = ts(kt, start = year[1], frequency = 1),  df = c(nu, nu.lin),
                   residuals = residuals, fitted = fitted, varprop = NA, # svd.mx$d[1]^2/sum(svd.mx$d^2),
                   y = fts(age, mx, start = year[1], frequency = 1, xname = "Age", yname = "Mortality"), 
                   mdev = c(sc.phi, sc.phi.lin), model = mod.choices[model], adjust = error, 
                   type = dat$type,
                   call = match.call(), conv.iter=iter, bx0=bx0, 
                   itx = ts(itx, start=year[1]-age[k], frequency=1))
    names(output)[4] <- if (is.numeric(series)) names(dat$rate)[series] else series
    names(output$mdev) <- c("Mean deviance base", "Mean deviance total")
    names(output$df) <- c("df base", "df tot")
    return(structure(output, class = c(paste("rh", model, sep=''), "rh", "lca", "fmm")))
}
