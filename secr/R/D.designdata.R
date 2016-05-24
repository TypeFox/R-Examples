###############################################################################
## package 'secr'
## D.design.MS.R
## Prepare density design matrix
##
## 2010 02 25 force levels=sessionlevels for dframe$session
## 2011 10 09 revised for multi-session mask covariates; needs testing
## 2011 11 08 revised meanSD, getcol, scale
## 2014-10-25 designmatrix
###############################################################################
## NOTE does not standardize sessioncov, maskcov
###############################################################################


###############################################################################
## used only for secondary setup (e.g. confint, score.test)

designmatrix <- function (modelled, mask, model, grouplevels, sessionlevels, 
                          sessioncov, smoothsetup) {
    if (!modelled) {
        designmatrix <- matrix(nrow = 0, ncol = 0)
        attr(designmatrix, 'dimD') <- NA
    }
    else {
        temp <- D.designdata( mask, model, grouplevels,
                              sessionlevels, sessioncov)
        designmatrix <- general.model.matrix(model, temp, smoothsetup)
        attr(designmatrix, 'dimD') <- attr(temp, 'dimD')
    }
    designmatrix
}
###############################################################################

D.designdata <- function (mask, Dmodel, grouplevels, sessionlevels, sessioncov = NULL,
                          meanSD = NULL) {

    ## mask -- mask object or list of masks of the same length as sessionlevels
    ## Dmodel -- formula that may be constant ~1 or include
    ## any of the 'automatic' terms c('g','x','y','x2','y2','xy','session',
    ## 'Session') or user-supplied mask-level covariates
    ## grouplevels -- group.levels(capthist,groups)
    ## sessionlevels -- character vector of session names
    ## sessioncov -- dataframe of session-specific covariates
    ## meanSD -- optional externally provided mean and SD (rows) for x-
    ##     and y-coordinates (columns)

    ## Output is a dataframe with one row for each combination of
    ## mask point, group and session. Conceptually, we use a rectangular
    ## array with enough rows to accommodate the largest mask, so some rows
    ## in the output may merely hold space to enable easy indexing.

    #--------------------------------------------------------------------------
    ## utility function

    ## standardise a numeric column from msk and pad with NA to desired length
    getcol <- function (msk, colnum) {
          mn <- attr(msk, "meanSD")[1, colnum]
          SD <- attr(msk, "meanSD")[2, colnum]
          pad1 (scale(msk[,colnum], mn, SD), nmaskrow)
    }

    #--------------------------------------------------------------------------
    ## setup
    vars  <- all.vars(Dmodel)
    ngrp  <- length(grouplevels)
    R     <- length(sessionlevels)
    if (ms(mask)) {
        maskrows <- sapply(mask, nrow)
        if (length(maskrows) != R)
            stop ("number of masks does not match number of sessions")
    }
    else {
        maskrows <- nrow(mask)
    }
    nmaskrow <- max(maskrows)
    dims  <- c(nmaskrow, ngrp, R)
    ## special case where new meanSD passed
    if (!is.null(meanSD)) {
        if (ms(mask))
            for (i in 1:length(mask))
                attr(mask[[i]], "meanSD") <- meanSD[[i]]
        else
            attr(mask, "meanSD") <- meanSD
    }

    #--------------------------------------------------------------------------
    ## coordinates
    ## might be condensed by operating on x,y together, but would it be clear?
    if (any (vars %in% c('x','y','x2','y2','xy'))) {
        if (ms(mask)) {
            ## session-specific masks
            x <- lapply(mask, getcol, 1)
            y <- lapply(mask, getcol, 2)
            x <- lapply(x, rep, ngrp)
            y <- lapply(y, rep, ngrp)
            dframe <- data.frame(
                x = as.vector(unlist(x)),
                y = as.vector(unlist(y))
            )
        }
        else {
            ## uniform mask across sessions
            x <- getcol(mask,1)
            y <- getcol(mask,2)
            dframe <- as.data.frame( list (
                x = rep(as.vector(unlist(x)), ngrp * R),
                y = rep(as.vector(unlist(y)), ngrp * R)
            ))
        }
        #---------------------------------------------
        ## coordinates transformed for quadratic trend
        if ('x2' %in% vars) {
            dframe$x2 <- dframe$x^2
        }
        if ('y2' %in% vars) {
            dframe$y2 <- dframe$y^2
        }
        if ('xy' %in% vars) {
            dframe$xy <- dframe$x * dframe$y
        }
    }
    else
        dframe <- data.frame(intercept = rep(1, nmaskrow * ngrp * R))

    #--------------------------------------------------------------------------
    ## groups
    if ('g' %in% vars) {
        if (length(grouplevels)<1)
            stop ("no groups specified")
        dframe$g <- insertdim(factor(grouplevels), 2, dims)
    }
    #--------------------------------------------------------------------------
    ## sessions
    if ('session' %in% vars) {
       dframe$session <- insertdim(factor(sessionlevels, levels =
           sessionlevels), 3, dims)
    }
    if ('Session' %in% vars) {
       dframe$Session <- insertdim(0:(R-1), 3, dims)
    }
    #--------------------------------------------------------------------------
    ## all autovars should have now been dealt with
    vars <- vars[!vars %in% c('g', 'x', 'y', 'x2', 'y2', 'xy',
        'session', 'Session')]

    #--------------------------------------------------------------------------
    ## session covariates
    if (!is.null(sessioncov)) {
        found <- names(sessioncov) %in% vars
        if (is.data.frame(sessioncov) & any(found)) {
            found <- names(sessioncov)[found]
            values <- as.data.frame(sessioncov[,found])
            names(values) <- found
            if (length(values)>0) {
                for (i in 1:ncol(values)) {
                    vals <- values[,i]
                    dframe[,found[i]] <- insertdim (vals, 3, dims)
                }
                vars <- vars[!(vars %in% found)]
            }
        }
    }
    #--------------------------------------------------------------------------
    ## handling of mask covariates substantially revised 2011-10-09
    maskcovs <- covariates(mask)
    if ((!is.null(maskcovs)) & (length(vars>0))) {
        std <- FALSE    ## no standardization
        ## requires all covariates in all masks
        expand <- function (df, n) {
            found <- names(df) %in% vars
            temp <- df[, found, drop = FALSE]
            ## pad with replicated row 1 to ensure rectangular
            if (nrow(temp) < n)
                rbind(temp, temp[rep(1,n-nrow(df)),,drop=FALSE])  ## drop=F added 2013-03-10
            else
                temp
        }
        if (ms(mask)) {
            maskcovs <- lapply(maskcovs, expand, nmaskrow)
            ncov <- sapply (maskcovs, ncol)
            if (any(ncov != ncov[1]))
                stop ("covariate missing in at least one session mask")
            maskcovs <- do.call(rbind, maskcovs)
        }
        else {
            maskcovs <- expand(maskcovs, nmaskrow)
        }
        for (i in names(maskcovs)) {
            vals <- maskcovs[,i] ## vector
            if (!is.factor(vals) & std)
                vals <- scale(vals)
            if (ms(mask))
                ## vals already full length to match nrow(dframe)
                dframe[,i] <- vals
            else
                ## vals repeated across groups & sessions
                dframe[,i] <- insertdim (vals, 1, dims)
        }
        vars <- vars[!(vars %in% names(maskcovs))]
    }
    #--------------------------------------------------------------------------
    ## unmatched variables in model?
    if (length(vars) > 0)
        stop (paste(vars,collapse=','), " not found")

    #--------------------------------------------------------------------------
    ## report dimensions as an attribute
    attr(dframe, 'dimD') <- c(nmaskrow, ngrp, R)
    attr(dframe, 'validMaskRows') <- maskrows
    dframe
}
###############################################################################
