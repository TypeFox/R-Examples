
".labDive" <- function(act, string)
{
    ## Value: Label dives along vector of same length as input.  Return a
    ## matrix labelling each dive and postdive reading
    ## --------------------------------------------------------------------
    ## Arguments: act=factor with values to label, string=character string
    ## to search in act to be labelled sequentially
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    dive <- vector(mode="numeric", length=length(act))
    dive[act == string] <- 1
    runs <- rle(dive)
    rawid <- rep(seq(along=runs$lengths), runs$lengths)

    diveid <- rawid
    diveid[dive == 0] <- 0             # non-dives are 0, and dives conseq:
    diveid[dive != 0] <- rep(seq(along=table(diveid)[-1]), table(diveid)[-1])
    pdid <- numeric(length(rawid))      # dives are 0 and postdives conseq:
    pdinds <- rawid %in% (unique(rawid[dive == 1]) + 1)
    pdid[pdinds] <- rep(seq(along=table(rawid[pdinds])), table(rawid[pdinds]))
    cbind(dive.id=diveid, postdive.id=pdid)
}


".detDive" <- function(zdepth, act, dive.thr)
{
    ## Value: A data frame; detecting dives, using a depth threshold
    ## --------------------------------------------------------------------
    ## Arguments: zdepth=depth vector of zoc'ed data, act=factor with
    ## dry/wet activity IDs (2nd element returned by detPhase), with values
    ## "W" for at-sea, dive.thr=dive threshold in m
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    ## Get the indices of below surface activity and label as "U"
    underw <- which((act == "W" | act == "Z") & zdepth > 0)
    if (length(underw) > 0) {
        act[underw] <- "U"
    } else warning("no subaquatic activity found")

    ## Label underwater excursions
    labuw <- .labDive(act, "U")
    ## Max depth of each "U" phase
    uwmax <- tapply(zdepth[underw], labuw[underw, 1], max, na.rm=TRUE)
    ## Change each "U" (underwater) phase to "D" (diving) if its max depth
    ## > dive threshold; i.e. we assume a dive started when depth dipped
    ## below zero and ended when it returned to zero *if the maximum depth
    ## of such a phase exceeded the dive threshold*.  This is a problem if
    ## we have noise above the dive threshold, so we need to switch the
    ## observations <= surface back to "U" and relabel the activity vector
    ## accordingly (this works great!)
    act[labuw[, 1] %in% as.numeric(names(uwmax[uwmax > dive.thr]))] <- "D"

    dives.maybe <- .labDive(act, "D")
    surface.idx <- which(dives.maybe[, 1] > 0 & zdepth <= dive.thr)
    act[surface.idx] <- "U"
    inddive <- .labDive(act, "D")
    ndives <- length(unique(inddive[act == "D", 1]))
    message(ndives, " dives detected")

    ## Return data frame with vectors of dive indices, adjusted activity,
    ## and postdive indices
    data.frame(dive.id=inddive[, 1], dive.activity=act,
               postdive.id=inddive[, 2])
}

##_+ Dive Detection with smoothing spline and derivative
".cutDive" <- function(x, smooth.par=NULL, knot.factor, descent.crit.q,
                       ascent.crit.q)
{
    ## Value: 'diveModel' object with details of dive phase model.
    ## --------------------------------------------------------------------
    ## Arguments: x=a 3-col numeric matrix with index in original TDR
    ## object, non-NA depths and numeric time.  A single dive's data
    ## (*below* 'dive.threshold'); smooth.par=spar parameter for
    ## smooth.spline(); knot.factor=numeric scalar that multiplies the
    ## duration of the dive (used to construct the time predictor for the
    ## derivative); descent.crit.q=ascent.crit.q quantiles defining the
    ## critical vertical rates of descent and ascent where descent should
    ## and ascent begin.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    ## We're passed sub-dive threshold observations, so we add 0 at ends to
    ## get proper derivatives
    depths <- c(0, x[, 2], 0)
    step.l <- if (nrow(x) < 2) 1 else median(diff(x[, 3]), na.rm=TRUE)
    times <- c(x[1, 3] - step.l, x[, 3], x[nrow(x), 3] + step.l)
    if (length(times) >= 4) {   # guard against smooth.spline() limitations
        times.scaled <- times - times[1]
    } else {
        times.scaledOrig <- times - times[1]
        times4 <- seq(times[1], times[length(times)], length.out=4)
        times4.scaled <- as.numeric(times4) - as.numeric(times4[1])
        depths.itpfun <- approxfun(times.scaledOrig, depths)
        depths <- depths.itpfun(times4.scaled)
        times.scaled <- seq(times.scaledOrig[1],
                            times.scaledOrig[length(times.scaledOrig)],
                            length.out=4)
    }
    ## Interpolate wiggles at ends; assume they are ZOC errors
    if ((length(times) > 4) && (depths[3] <= depths[2])) {
        depths[2] <- approx(1:2, c(depths[3], 0), xout=1.5)$y
    }
    if ((length(times) > 4) &&
        (depths[length(depths) - 1] < max(depths)) &&
        (depths[length(depths) - 1] >= depths[length(depths) - 2])) {
        depths[length(depths) - 1] <- approx(1:2,
                                             c(depths[length(depths) - 2], 0),
                                             xout=1.5)$y
    }
    times.pred <- seq(times.scaled[1], times.scaled[length(times.scaled)],
                      length.out=length(times.scaled) * knot.factor)
    if (is.null(smooth.par)) {
        spl <- stats::smooth.spline(times.scaled, depths, all.knots=TRUE)
        depths.smooth <- stats::smooth.spline(times.scaled, depths,
                                              spar=spl$spar, all.knots=TRUE)
    } else {
        depths.smooth <- stats::smooth.spline(times.scaled, depths,
                                              spar=smooth.par, all.knots=TRUE)
    }
    depths.deriv <- predict(depths.smooth, times.pred, deriv=1)
    depths.d <- depths.deriv$y

    ## Descent ------------------------------------------------------------
    Dd1pos.rle <- rle(sign(depths.d))
    ## First sequence of positives (descent)
    Dd1pos.idx <- which.max(Dd1pos.rle$values)
    ## How long is this sequence from beginning?
    Dd1pos.sum <- sum(Dd1pos.rle$lengths[seq(Dd1pos.idx)])
    Dd1.maybe <- depths.d[seq(Dd1pos.sum)]
    ## Index with first minimum from beginning (only positives)
    if (all(Dd1.maybe <= 0)) {     # but first maximum if all non-positives
        Dd1pos.min <- which.max(Dd1.maybe)
    } else {
        d.crit.rate <- quantile(Dd1.maybe, probs=descent.crit.q)
        beyond <- Dd1.maybe > d.crit.rate
        Dd1pos.min <- ifelse(any(beyond),
                             which.min(Dd1.maybe[beyond]),
                             which.min(Dd1.maybe >= 0))
    }
    ## Or one could just do this if not worried about starting negative d'
    ## Dd1pos.l <- rle(sign(depths.d))$lengths[1]
    ## Dd1pos.min <- which.min(depths.d[seq(Dd1pos.l)])
    Dd1pos.crit <- which.min(abs(times.pred[Dd1pos.min] - times.scaled))

    ## Ascent -------------------------------------------------------------
    Ad1 <- rev(depths.d)                # work from end of dive
    Ad1neg.rle <- rle(sign(Ad1))
    ## First sequence of negatives (ascent)
    Ad1neg.idx <- which.min(Ad1neg.rle$values)
    ## How long is this sequence from the end?
    Ad1neg.sum <- sum(Ad1neg.rle$lengths[seq(Ad1neg.idx)])
    ## Potential ascent derivatives from end
    Ad1.maybe <- Ad1[seq(Ad1neg.sum)]
    ## Potential ascent derivatives in natural order
    Ad1.maybe.nat <- rev(Ad1[seq(Ad1neg.sum)])
    ## Index with first maximum in natural ascent order (only negatives)
    if (all(Ad1.maybe >= 0)) {      # but first minimum if all non-negative
        Ad1neg.max.nat <- which.min(Ad1.maybe.nat)
    } else {
        a.crit.rate <- quantile(Ad1.maybe.nat, probs=(1 - ascent.crit.q))
        beyond <- Ad1.maybe.nat < a.crit.rate
        beyond.w <- which(beyond)    # indices below critical in candidates
        beyond0 <- Ad1.maybe.nat < 0
        beyond0.w <- which(beyond0) # indices of non-positives in candidates
        Ad1neg.max.nat <- ifelse(any(beyond),
                                 beyond.w[which.max(Ad1.maybe.nat[beyond])],
                                 beyond0.w[which.max(Ad1.maybe.nat[beyond0])])
    }
    ## Position of *first* maximum negative derivative in the sequence at
    ## the end of the predicted smooth depths (from end)
    Ad1neg.max <- Ad1neg.sum - Ad1neg.max.nat + 1
    ## Absolute differences between time predictor corresponding to the
    ## position above and scaled time (it's important it's in reverse, so
    ## that we can find the first minimum below)
    Ad1neg.diff <- rev(times.pred)[Ad1neg.max] - rev(times.scaled)
    left <- max(which(Ad1neg.diff <= 0))
    right <- left + 1
    left.idx <- length(times.scaled) - left + 1
    right.idx <- length(times.scaled) - right + 1
    both.idx <- length(times.scaled) - which.min(abs(Ad1neg.diff)) + 1
    ## Position of the *first* minimum absolute difference above
    Ad1neg.crit <- ifelse(depths[left.idx] == depths[right.idx],
                          left.idx, both.idx)

    ## Correct both critical indices if this is a brief dive
    if (length(times) < 4) {
        Dd1pos.crit <- which.min(abs(times.scaled[Dd1pos.crit] -
                                     times.scaledOrig))
        Ad1neg.crit <- which.min(abs(times.scaled[Ad1neg.crit] -
                                     times.scaledOrig))
    }
    ## Correct for added 0 at ends
    descind <- seq(Dd1pos.crit - 1)
    ascind <- seq(min(nrow(x), Ad1neg.crit - 1), length(times) - 2)

    ## Bottom -------------------------------------------------------------
    bottind <- c(descind[length(descind)],
                 setdiff(seq(nrow(x)), union(descind, ascind)),
                 ascind[1])

    ## descent is everything in descind that's not in union of bottind and ascind
    d <- setdiff(descind, union(bottind, ascind))
    ## descent/bottom is what's common to descind and bottind
    db <- intersect(descind, bottind)
    ## bottom is everything in bottind that's not in union of descind and ascind
    b <- setdiff(bottind, union(descind, ascind))
    ## bottom/ascent is what's common to ascind and bottind
    ba <- intersect(ascind, bottind)
    ## descent/ascent is what's common to descind and ascind
    da <- intersect(descind, ascind)
    ## ascent is everything in ascind that's not in union of descind and bottind
    a <- setdiff(ascind, union(descind, bottind))

    labs <- character(nrow(x))
    labs[d] <- "D"
    labs[db] <- "DB"
    labs[b] <- "B"
    labs[ba] <- "BA"
    labs[da] <- "DA"
    labs[a] <- "A"
    rowids <- unique(c(x[d, 1], x[db, 1], x[b, 1], x[ba, 1], x[a, 1], x[da, 1]))
    ## Any remaining rowids are nonexistent labels
    label.mat <- cbind(rowids[!is.na(rowids)], labs[!is.na(rowids)])
    new("diveModel",
        label.matrix=label.mat,
        dive.spline=depths.smooth,
        spline.deriv=depths.deriv,
        descent.crit=Dd1pos.crit,
        ascent.crit=Ad1neg.crit,
        descent.crit.rate=d.crit.rate,
        ascent.crit.rate=a.crit.rate)
}

".labDivePhase" <- function(x, diveID, ...)
{
    ## Value: A list with a factor labelling portions of dives, and a list
    ## of 'diveModel' objects for each dive.
    ## --------------------------------------------------------------------
    ## Arguments: x=class TDR object, diveID=numeric vector indexing each
    ## dive (non-dives should be 0). As it is called by calibrateDepth,
    ## these indices include underwater phases, not necessarily below dive
    ## threshold. ...=arguments passed to .cutDive(), usually from
    ## calibrateDepth() and include 'smooth.par' and 'knot.factor'.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    if (!is(x, "TDR")) stop("x must be a TDR object")
    ok <- which(diveID > 0 & !is.na(getDepth(x))) # required diving indices

    if (length(ok) > 0) {
        ddepths <- getDepth(x)[ok]               # diving depths
        dtimes <- getTime(x)[ok]                 # diving times
        dids <- diveID[ok]                       # dive IDs
        ## We send a matrix of indices, and non-NA depths and times
        td <- matrix(data=c(ok, ddepths, as.numeric(dtimes)), ncol=3)
        ## Problems with by() always returning 'by', not list
        perdivetd <- c(by(td, dids, .cutDive, ...,
                          simplify=FALSE))
        labdF <- do.call(rbind, lapply(perdivetd, slot, "label.matrix"))
        ff <- factor(rep("X", length(diveID)),
                     levels=c(unique(labdF[, 2]), "X"))
        ff[as.numeric(labdF[, 1])] <- labdF[, 2]
        list(phase.labels=ff, dive.models=perdivetd)
    } else {
        warning("no dives were found in x, hence no dive models created")
        list(phase.labels=factor(rep("X", length(diveID))),
             dive.models=list())
    }
}

".diveMatches" <- function(diveID, diveNo)
{
    ## Value: Logical vector with elements of diveNo having a match in
    ## diveID.
    ## --------------------------------------------------------------------
    ## Arguments: diveID=numeric or character vector with purported dive
    ## numbers; diveNo=Numeric vector with dive numbers that do exist in
    ## diveID.
    ## --------------------------------------------------------------------
    ## Purpose: Check for the existence of dive numbers (diveNo) in a list
    ## (l).
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    dinl <- diveNo %in% diveID
    if (! any(dinl)) stop("none of the requested dives exist")
    if (any(! dinl))
        warning("dives ", diveNo[! dinl],
                " do(es) not exist, so is(are) ignored")
    dinl
}

".diveIndices" <- function(diveID, diveNo)
{
    ## Value: A numeric vector with the indices of dives (and their
    ## beginning/end indices) in diveID
    ## --------------------------------------------------------------------
    ## Arguments: diveID=numeric vector numbering all dives and non-dives,
    ## diveNo=numeric vector of unique dive indices to extract fromdiveID.
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    tryCatch({
        didx <- .diveMatches(diveID, diveNo)
        ok <- which(diveID %in% diveNo[didx])
        okl <- pmax(1, setdiff(ok - 1, ok))
        okr <- pmin(length(diveID), setdiff(ok + 1, ok))
        unique(sort(c(okl, ok, okr)))   # add the surface points
    })
}


## TEST ZONE --------------------------------------------------------------

## utils::example("calibrateDepth", package="diveMove", ask=FALSE, echo=FALSE)
## X <- c(2, 7, 100, 120, 240)
## ## .labDivePhase(getTDR(tdr.calib), getDAct(tdr.calib, "dive.id"),
## ##                          smooth.par=0.1, knot.factor=30, descent.crit=0.01,
## ##                          ascent.crit=0)
## ## diveX <- as.data.frame(extractDive(dcalib, diveNo=X[5]))
## X <- c(2, 7, 100, 120, 743, 1224, 1222, 1223)
## diveX <- as.data.frame(extractDive(tdr.calib, diveNo=X[8]))
## diveX.m <- cbind(as.numeric(row.names(diveX[-c(1, nrow(diveX)), ])),
##                  diveX$depth[-c(1, nrow(diveX))],
##                  diveX$time[-c(1, nrow(diveX))])
## phases <- .cutDive(diveX.m, smooth.par=0.1, knot.factor=30,
##                    descent.crit.q=0.01, ascent.crit.q=0)

## for (dive in seq(max(dives[dives > 0]))) {
##     diveX <- as.data.frame(extractDive(tdr.calib, diveNo=dive))
##     diveX.m <- cbind(as.numeric(row.names(diveX[-c(1, nrow(diveX)), ])),
##                      diveX$depth[-c(1, nrow(diveX))],
##                      diveX$time[-c(1, nrow(diveX))])
##     phases <- .cutDive(diveX.m, smooth.par=0.1, knot.factor=50,
##                        descent.crit.q=0.01, ascent.crit.q=0)
## }
