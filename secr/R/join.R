## 2011-12-01, 2011-12-07, 2012-06-08, 2012-11-04
## join returns single-session object from list of inputs

## should new trap ID be numeric , character or factor???
## 2015-01-31 join() bug fix could fail with exclusivedetector types
## 2015-04-05 previous bug fix failed...
## 2015-04-13
## 2015-10-29 extended to deal with Tu, Tm and all-zero histories
## 2016-01-08 addzerodf function moved to utility.R

join <- function (object, remove.dupl.sites = TRUE, tol = 0.001) {

    ####################################################################
    onesession <- function (sess) {
        ## form CH as a dataframe
        CH <- object[[sess]]
        newID <- animalID(CH)
        newocc <- occasion(CH) + before[sess]
        newtrap <- trap(CH)
        df <- data.frame(newID = newID, newocc = newocc, newtrap = newtrap,
                         alive = alive(CH), sess = rep(sess, length(newID)),
                         stringsAsFactors = FALSE)
        if (!is.null(xy(CH)))
            df[,c('x','y')] <- xy(CH)
        if (!is.null(signal(CH)))
            df[,'signal'] <- signal(CH)
        
        ## add all-zero (sighting-only) records as required,
        ## otherwise leave unchanged
        addzerodf (df, CH, sess)
    }
    ####################################################################

    condition.usage <- function (trp, i) {
        us <- matrix(0, nrow=nrow(trp), ncol=nnewocc)
        s1 <- c(1, cumsum(nocc)+1)[i]
        s2 <- cumsum(nocc)[i]
        if (is.null(usage(trp)))
            us[,s1:s2] <- 1
        else
            us[,s1:s2] <- usage(trp)
        usage(trp) <- us
        trp
    }
    ####################################################################
    ## preparing for merge when traps vary... 2015-10-29
    condition.sightings <- function (CH, i, type = 'Tu') {
        T <- attr(CH, type)
        if (!is.null(T)) {
            if (is.matrix(T)) {
                Tnew <- matrix(0, nrow=nrow(traps(CH)), ncol=nnewocc)
                s1 <- c(1, cumsum(nocc)+1)[i]
                s2 <- cumsum(nocc)[i]
                Tnew[,s1:s2] <- T
                attr(CH, type) <- Tnew
            }
        }
        CH
    }
    ####################################################################

    ## mainline

    if (!ms(object) | any(sapply(object, class) != 'capthist'))
        stop("requires multi-session capthist object or list of ",
             "single-session capthist")

    outputdetector <- detector(traps(object)[[1]])
    nsession <- length(object)
    nocc <- sapply(object, ncol)
    names(nocc) <- NULL
    nnewocc <- sum(nocc)
    ## cumulative number of preceding occasions
    before <- c(0, cumsum(nocc)[-nsession])

    ##------------------------------------------------------------------
    ## combine capthist as one long dataframe
    df <- lapply(1:nsession, onesession)
    df <- do.call(rbind, df)
    n <- length(unique(df$newID))

    ##------------------------------------------------------------------
    ## resolve traps
    ## first check whether all the same (except usage)
    temptrp <- lapply(traps(object), function(x) {usage(x) <- NULL; x})
    sametrp <- all(sapply(temptrp[-1], identical, temptrp[[1]]))

    if (sametrp & remove.dupl.sites) {
        newtraps <- temptrp[[1]]
        class(newtraps) <- c("traps", "data.frame")
        if (length(usage(traps(object))) > 0)
            usage(newtraps) <- do.call(cbind, usage(traps(object)))
        ## df$newtrap unchanged
    }
    else {
        ## temptrp <- mapply(condition.usage, traps(object), 1:nsession)
        ## 2015-12-08
        temptrp <- mapply(condition.usage, traps(object), 1:nsession, SIMPLIFY = FALSE)
        newtraps <- do.call(rbind, c(temptrp, renumber = FALSE))
        class(newtraps) <- c("traps", "data.frame")
        df$newtrap <- paste(df$newtrap,df$sess, sep=".")
    }

    ##------------------------------------------------------------------
    ## ensure retain all occasions
    df$newocc <- factor(df$newocc, levels = 1:nnewocc)

    ##------------------------------------------------------------------
    ## construct new capthist matrix or array from positive detections
    if (outputdetector %in% .localstuff$exclusivedetectors) {
        alivesign <- df$alive*2 - 1
        tempnew <- matrix(0, nrow = n, ncol = sum(nocc))
        dimnames(tempnew) <- list(unique(df$newID), 1:sum(nocc))

        ## bug fix 2015-04-05
        ## trapindex <- match(df$newtrap, unique(df$newtrap))
        trapindex <- match(df$newtrap, rownames(newtraps))
        tempnew[cbind(df$newID, df$newocc)] <- trapindex * alivesign
    }
    else {
        if (outputdetector %in% .localstuff$polydetectors)
            df$newtrap <- factor(df$newtrap)
        else
            df$newtrap <- factor(df$newtrap, levels=rownames(newtraps))
        tempnew <- table(df$newID, df$newocc, df$newtrap, useNA = "no")
        alivesign <- tapply(df$alive, list(df$newID,df$newocc,df$newtrap),all)
        alivesign[is.na(alivesign)] <- TRUE
        alivesign <- alivesign * 2 - 1
        tempnew <- tempnew * alivesign
    }

    ##------------------------------------------------------------------
    ## pile on the attributes...
    class(tempnew) <- 'capthist'
    traps(tempnew) <- newtraps
    session(tempnew) <- 1
    neworder <- order (df$newocc, df$newID, df$newtrap)

    ##------------------------------------------------------------------
    ## concatenate marking-and-resighting-occasion vectors
    tempmarkocc <- unlist(markocc(traps(object)))
    if (!is.null(tempmarkocc)) {
        names(tempmarkocc) <- NULL
        markocc(traps(tempnew)) <- tempmarkocc
    }
    
    ##------------------------------------------------------------------
    ## unmarked and nonID sightings
    ## not yet implemented for varying traps
    if (sametrp & remove.dupl.sites) {
        ## retain unmarked sightings and nonID sightings if present
        ## ignore if NULL
        Tu <- Tu(object)
        if (!is.null(Tu[[1]])) {
            if (!all(sapply(Tu, is.matrix)))
                Tu(tempnew) <- do.call(sum, Tu)
            else
                Tu(tempnew) <- do.call(cbind, Tu)
        }

        Tm <- Tm(object)
        if (!is.null(Tm[[1]])) {
            if (!all(sapply(Tm, is.matrix)))
                Tm(tempnew) <- do.call(sum, Tm)
            else
                Tm(tempnew) <- do.call(cbind, Tm)
        }
    }
    else {
        ## Tu, Tm not ready yet
        if (!is.null(Tu(object[[1]])) | !is.null(Tm(object[[1]])))
            stop ("join does not yet merge sighting matrices when traps vary")
    }

    ##------------------------------------------------------------------
    ## covariates, xy, signal attributes
    if (!is.null(covariates(object))) {
        tempcov <- do.call(rbind, covariates(object))
        if (!is.null(tempcov)) {
            IDcov <- unlist(lapply(object,rownames))
            ## use first match
            tempcov <- tempcov[match(rownames(tempnew), IDcov),,drop = FALSE]
            rownames(tempcov) <- rownames(tempnew)
            covariates(tempnew) <- tempcov
        }
    }

    ##------------------------------------------------------------------
    ## negotiate problem that all-zero histories have no xy, signal
    tempdf <- df[neworder,, drop = FALSE]
    if (!is.null(df$x)) {
        xy(tempnew) <- tempdf[!is.na(tempdf$newocc), c('x','y')]
    }
    if (!is.null(df$signal))
        signal(tempnew) <- tempdf[!is.na(tempdf$newocc),'signal']

    ##------------------------------------------------------------------
    ## purge duplicate sites, if requested
    if (remove.dupl.sites & !sametrp)
        tempnew <- reduce(tempnew, span=tol, dropunused = FALSE, verify = FALSE)

    ## remember previous structure, for MARK-style robust design
    attr(tempnew, 'interval') <- unlist(sapply(nocc, function(x) c(1,rep(0,x-1))))[-1]

    ##------------------------------------------------------------------

    tempnew

}

RMarkInput <- function (object, grouped = FALSE, covariates = TRUE) {
    if (!inherits(object, "capthist"))
        stop ("requires single-session capthist object")
    if (ms(object))
        stop ("requires single-session capthist object - use 'join'")
    if (length(dim(object)) != 2)
        CH0 <- reduce(object, outputdetector = 'multi')
    else
        CH0 <- object
    ntimes <- ncol(object)
    alive <- apply(object,1,function(x) all(x>=0))
    CH <- pmin(abs(CH0),1)

    if (is.logical(covariates)) {
        if (covariates) {
            if (is.null(covariates(CH0)))
                stop("no covariates in object")
            covnames <- names(covariates(CH0))
        }
        else
            covnames <- ""
    }
    else {
        covnames <- covariates
        if (is.character(covariates)) {
            if (is.null(covariates(CH0)))
                stop("no covariates in object")
        }
        found <- covnames %in% names(covariates(CH0))
        if (any(!found)) {
            stop(paste(covnames[!found], collapse=','), " not in covariates(object)")
        }
    }

    if (any(covnames != "")) {
        if (grouped) {
            warning("'grouped' is incompatible with individual covariates and will be ignored")
            grouped <- FALSE
        }
    }

    if (grouped)   ## bug fix 2012-07-04
        CH <- cbind(CH, alive) ## add single-digit code as last column
    CH <- data.frame(ch = apply(CH, 1, paste, collapse=''),
        stringsAsFactors = FALSE)

    if (grouped) {
        temp <- table(CH$ch)
        alive <- as.numeric(substring(names(temp),ntimes+1,ntimes+1))
        CH <- data.frame(ch = substring(names(temp),1,ntimes),
                         freq = as.numeric(temp),
                         stringsAsFactors = FALSE)
        CH$freq <- ifelse(alive, CH$freq, -CH$freq)
        CH <- CH[order(CH$ch, decreasing = TRUE),]
        row.names(CH) <- 1:nrow(CH)
    }
    else {
        CH$freq <- ifelse(alive,1,-1)
        if (any(covnames != "")) {
            CH[,covnames] <- covariates(CH0)[,covnames]
        }
    }
    attr(CH, "interval") <- attr(object, "interval")
    if (is.null(attr(CH,"interval")))
        attr(CH, "interval") <- rep(0,ntimes-1)
    CH
}

# temp <- join(ovenCH)
# ovenint <- c(rep(0,8),rep(c(1,rep(0,9)),4))
# secr.fit(temp, details=list(interval = ovenint))

# plot(join(ovenCH), tracks=T)

unjoin <- function (object, interval, ...) {
    if (missing(interval) & is.null(attr(object,'interval')))
        stop ("requires 'interval' to define sessions")
    if (missing(interval) )
        interval <- attr(object,"interval")
    session <- c(0,cumsum(interval>0))+1
    nsess <- max(session)
    if (nsess<2) {
        warning ("interval define only one session")
        return(object)
    }
    newobj <- vector(mode='list', length=nsess)
    for (sess in 1:nsess) {
        newobj[[sess]] <- subset(object, occasions = (session==sess), ...)
    }
    class (newobj) <- c('list','capthist')
    session(newobj) <- 1:nsess
    return(newobj)
}

unRMarkInput <- function(df, covariates = TRUE) {
    if (!is.data.frame(df))
        stop("requires dataframe input")
    if (!all(c('ch','freq') %in% names(df)))
        stop ("ch and freq are required fields")
    nocc <- nchar(df$ch)
    if (length(unique(nocc))>1)
        stop ("ch must be a constant-length string of 0s and 1s")
    nocc <- nocc[1]
    freq <- df$freq
    alive <- sign(freq)
    freq <- abs(freq)
    freq <- rep(1:nrow(df), freq)
    alive <- alive[freq]
    ch <- df$ch[freq]
    CH <- matrix(as.numeric(unlist(sapply(ch, strsplit, ''))), byrow = TRUE, ncol = nocc)
    class(CH) <- 'capthist'
    # allow deads
    last <- function(x) which.max(cumsum(x))
    CH[cbind(1:nrow(CH), apply(CH,1,last))] <- alive
    # transfer covariates if present
    if (ncol(df)>2) {
        if (is.logical(covariates)) {
            if (covariates)
                covnames <- names(df)[-match(c('ch','freq'),names(df))]
            else
                covnames <- ""
        }
        else {
            covnames <- covariates[covariates %in% names(df)]
        }

        if (any(covnames != ""))
        covariates(CH) <- df[freq, covnames, drop = FALSE]
    }
    CH
}
