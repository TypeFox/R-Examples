elca.rh <-
function(dat,  year=dat$year, age=dat$age, dec.conv = 6,  
                    error = c("poisson", "gaussian"), 
                    restype = c("logrates", "rates", "deaths", "deviance"), scale = F, 
                    interpolate = F, verbose = T, spar=NULL, ax.fix = NULL){
    
    # --- Preamble: --- #
    # check data type:
    if (class(dat) != "rhdata") stop("Not \"rhdata\" class mortality data object!")
    
    # Check data compatibility:
    if (length(dat$covariates) > 1) stop('This function can only fit one extra paramter')
    
    cat('Original sample: ')
    print(dat)
    
    # Check the required age and period data ranges:
    ind <- age%in%dat$age
    if (!all(ind)) warning(paste('The following ages are not available:',
                                 paste(age[!ind], collapse=' '), sep='\n'))
    ind <- year%in%dat$year
    if (!all(ind)) warning(paste('The following years are not available:',
                                 paste(year[!ind], collapse=' '), sep='\n'))
    
    # Restrict data set to the required age and period.
    ind <- dat$age%in%age
    if (!all(ind)){
        dat$deaths <- dat$deaths[ind,,]
        dat$pop <- dat$pop[ind,,]
        dat$mu <- dat$mu[ind,,]
        dat$age <- dat$age[ind]
    }
    age <- dat$age
    ind <- dat$year%in%year
    if (!all(ind)){
        dat$deaths <- dat$deaths[,ind,]
        dat$pop <- dat$pop[,ind,]
        dat$mu <- dat$mu[,ind,]
        dat$year <- dat$year[ind]
    }
    year <- dat$year
    
    # Check error type:
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
    
    # first extract and check the force of mortality:
    mx <- dat$mu
    ind <- bool(mx < 1e-09, na=T) # get 0 (or negative) and NA value index
    if (sum(ind)) { 
        if (interpolate) {  
            # re-estimate all 0 and NA mu-s:
            dat <- fill.rhdata(dat, method='interp')
            mx <- dat$mu
            ind <- bool(mx < 1e-09, na=T) # reset index for final check
            if (sum(ind)) stop('It seems that interpolation did not help!\n')
        } else {
            # re-set all 0 mu-s to NA (due to log(0)=Inf )
            mx[ind] <- NA
            warning(paste("There are", sum(ind), "cells with 0/NA mu, which are", 
                          "ignored in the current analysis.\n  Try reducing the age range", 
                          "or setting interpolate=TRUE."))
        }
    }    
    mz <- log(apply(mx, 1:2, mean, na.rm=T))
    
    # Once the mu-s are checked, extract the exposure matrix and the weights:
    me <- dat$pop
    mw <- if (error == 'poisson') bool(me>1e-9, na=F) else array(1, dim=dim(mx)) 
    if (sum(!mw)){
        warning(paste("There are", sum(!mw), "cells with 0/NA exposures,", 
                      "which are ignored in the current analysis.\n  Try reducing the", 
                      "fitted age range.\n  Alternatively, fit ELC", 
                      "model with error=", mark('gaussian', F), "."))
    }
    
    cat('Applied sample: ')
    print(dat)
    
    # Finally, get the number of deaths
    md <- dat$deaths
    n <- length(year); k <- length(age)  # ; rtx <- n+k-1  # cohort range
    g <- length(dat$covariates[[1]])
    
    ft <- gl(n, k) 
    
    mod.choices <- c(' LC(g) = a(x)+a(g)+b(x)*k(t) ')    
    cat('\n  Fitting model:', sqb(mod.choices[1]), '\n\t- with', error, 
        'error structure', if (error=='poisson') 'and with deaths as weights', 
        '-\n')
    
    if (verbose){
        cat('Note:', sum(bool(md < 1e-9, na=T)), 'cells have 0/NA deaths and ',
            sum(!mw), 'have 0/NA exposure', '\n  out of a total of', n*k*g, 
            'data cells.\n')
    }
    
    # --- End Preamble. Start fitting: --- #
    
    # get starting values (unless ax is fixed):
    # this gives the log of the geometrical means of the mortality rates:
    ax <- if (is.null(ax.fix)) apply(mz, 1, mean, na.rm=T) else ax.fix
    ag <- rep(0, g); names(ag) <- dat$covariates[[1]]
    
    # fitting the model as a single vector using the defined factors (see setup)
    kt <- rep(0, n); names(kt) <- year
    bx <- rep(1/k, k); names(bx) <- age
    
    if (verbose) {
        cat('\n Starting values are:\n')
        print(structure(list(age=ax, int=bx, per=kt, ag), class='coef'))
    }
    
    # assign target (logrates or number of deaths):
    if (error == 'gaussian') { my <- mz } 
    if (error == 'poisson') { my <- md } 
    # fit LCM (Gaussian error structure, central rates as target)
    mf <- elc(ax, ag, bx, kt)
    # fit LCM (Poisson error structure, deaths as target)
    if (error == 'poisson') mf <- me*exp(mf)
    # Calculate deviance:
    D <- deviance.lca(my, mf, mw, error)    
    if (verbose) cat('\n Iterative fit:\n #iter   Dev    non-conv\n')
    dD <- 1 # a fictive difference of 1
    ncv <- iter <- 0  # non-convergence and iteration counters
    while(round(dD, dec.conv)){   # && iter <  end.iter
        # Note: with higher precisions it might over-fit the parameters 
        iter <- iter + 1
        if (verbose) cat('  ', c(iter, round(D, 6), ncv), '\n ', sep='  ')
        # --------  Stage 0: update ax ----------
        if (is.null(ax.fix)){
            td <- apply(mw*(my-mf), 1, sum, na.rm=T)
            edt <- mw
            if (error == 'poisson') edt <- edt*mf
            edt <- apply(edt, 1, sum, na.rm=T)
            ax <- ax + td/edt
            # fit LC (Gaussian error structure, central rates as target)
            mf <- elc(ax, ag, bx, kt)
            # fit LC (Poisson error structure, deaths as target)
            if (error == 'poisson') mf <- me*exp(mf)
        }
        # --------  Stage 1: update ag ----------
        td <- apply(mw*(my-mf), 3, sum, na.rm=T)
        edt <- mw
        if (error == 'poisson') edt <- edt*mf
        edt <- apply(edt, 3, sum, na.rm=T)
        ag <- ag + td/edt
        # adjust updates (to make them comparable to the glm coefficients):
        ag <-  ag-ag[1] 
        # fit LC (Gaussian error structure, central rates as target)
        mf <- elc(ax, ag, bx, kt)
        # fit LC (Poisson error structure, deaths as target)
        if (error == 'poisson') mf <- me*exp(mf)
        # --------  Stage 2: update kt ----------
        td <- apply(mw*(my-mf)*bx, 2, sum, na.rm=T)
        edt <- mw*bx^2
        if (error == 'poisson') edt <- edt*mf
        edt <- apply(edt, 2, sum, na.rm=T)
        kt <- kt + td/edt; 
        # apply sum(kt)=0 constrain to the period effect kt:
        kt <-  kt-mean(kt) 
        # fit LCM (Gaussian error structure, central rates as target)
        mf <- elc(ax, ag, bx, kt)
        # fit LCM (Poisson error structure, deaths as target)
        if (error == 'poisson') mf <- me*exp(mf)
        # --------  Stage 3: update bx ----------
        td <- apply(mw*(my-mf)*kt[ft], 1, sum, na.rm=T)
        edt <- mw*kt[ft]^2
        if (error == 'poisson') edt <- edt*mf
        edt <- apply(edt,  1, sum, na.rm=T)
        bx <- bx + td/edt
        if (!is.null(spar)) bx <- smooth.spline(bx, spar=spar)$y
        # fit LCM (Gaussian error structure, central rates as target)
        mf <- elc(ax, ag, bx, kt)
        # fit LCM (Poisson error structure, deaths as target)
        if (error == 'poisson') mf <- me*exp(mf)
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
    
    # reset the names of bx if smoothing applied:
    if (!is.null(spar)) names(bx) <- age
    
    cat('\n Iterations finished in:', iter, 'steps\n')
    # apply sum(bx)=1 constraint to the interaction parameter bx 
    sb <- sum(bx)
    bx <- bx/sb; kt <- kt*sb
    
    # Prepare return object:
    mf <- elc(ax, ag, bx, kt)  # fitted logrates
    mdf <- me*exp(mf)     # fitted number of deaths
    # for deviance residuals:
    mD <- deviance.lca(my, if (error=='poisson') mdf else mf, mw, error, total=F)
    # degrees of freedom depending on model and adjusted for missing data cells:
    nu <- (k-1)*(g-1)*(n-2) - sum(!mw)   
    sc.phi <- sum(mD, na.rm=T)/nu  # scaling factor (phi)
    # Note: the factor of 2 (given in the Hyndman function) is already included 
    # in the deviance.lca function.
    fit <- switch(restype, mf, exp(mf), mdf, if (error=='poisson') mdf else mf)
    fitted <- apply(fit, 3, fts, x=age, f = 1, start = year[1], xname = "Age", 
                    yname = paste("Fitted", rb(res.choices[if (restype<4) restype else 
                        if (error=='poisson') 3 else 1]), "mortality rate"))
    fitted$type <- res.choices[if (restype<4) restype else if (error=='poisson') 3 else 1]
    res <- switch(restype, log(mx), mx, md, my) - fit
    if (restype == 4) res <- sign(res)*sqrt(mD/sc.phi)    
    residuals <- apply(res, 3, fts, x=age, f = 1, start = year[1], xname = "Age",  
                       yname = paste("Residuals", rb(res.choices[restype]), "mortality rate"))
    residuals$type <- res.choices[restype]
    if (scale) {
        avdiffk <- -mean(diff(kt))
        bx <- bx * avdiffk
        kt <- kt/avdiffk
        warning('\"kt\" and \"bx\" parameters were re-scaled!')
    }
    if (verbose) {
        cat('  Updated values are:\n')
        print(structure(list(age=ax, int=bx, per=kt, ag), class='coef'), dec=5)
        # print out total sums too:
        cat('\t total sums are: \n')                
        psu <- list(bx=bx, kt=kt)
        psu <- unlist(lapply(psu, sum))
        print(round(psu, 4))
    }
    
    # Alternative deviance scaling by Hyndman (needs checking)
    drift <- mean(diff(kt))
    kt.lin <- mean(kt) + drift * (1:n - (n + 1)/2)
    mdf.lin <- me * exp(elc(ax, ag, bx, kt.lin))
    nu.lin <- k * g *(n - 2)   # for simple LC k/(k-1)*nu
    sc.phi.lin <- deviance.lca(md, mdf.lin, mw,'poisson')/nu.lin
    kt <- ts(kt, start = year[1], frequency = 1)
    y <- apply(mx, 3, fts, x=age, start = year[1], f = 1, xname = "Age", 
               yname = "Mortality")
    output <- list(); 
    output$lca <- vector('list', g)
    names(output$lca) <- dat$covariates[[1]]
    label <- paste(dat$label, dat$name, sep=' : ')
    labelX <- paste('X', dat$covariates[[1]], sep='.')
    for(i in seq(g)){
        output$lca[[i]] <- list(label = label, 
                                age = age, year = year, mx = mx[,,i], 
                                ax = ax+ag[i], bx = bx, kt = kt, adjust = error,
                                type = dat$type,
                                residuals = residuals[[i]], fitted = fitted[[i]],  y = y[[i]])
        names(output$lca[[i]])[4] <- labelX[i]
        output$lca[[i]] <- structure(output$lca[[i]], class = c("lca", "fmm")) 
    }
    output$age <- age
    output$year <- year
    output$ag <- ag
    output$ax <- ax
    output$bx <- bx   
    output$kt <- kt
    output$adjust <- error
    output$type <- dat$type
    output$label <- label
    output$call <- match.call()
    output$conv.iter <- iter
    output$varprop <- NA # svd.mx$d[1]^2/sum(svd.mx$d^2),
    output$mdev <- c(sc.phi, sc.phi.lin)
    names(output$mdev) <- c("Mean deviance base", "Mean deviance total")
    output$model <- mod.choices[1]
    output$df <- c(nu, nu.lin)
    names(output$df) <- c("df base", "df tot")
    return(structure(output, class=c("elca")))
}
