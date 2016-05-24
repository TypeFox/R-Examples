
###_ + Show and plot

###_  . TDR and TDRspeed
setMethod("show", signature=signature(object="TDR"),
          definition=function(object) {
              trange <- range(object@time)
              cat("Time-Depth Recorder data -- Class",
                  class(object), "object\n")
              cat("  Source File          : ", object@file, "\n",
                  sep="")
              cat("  Sampling Interval (s): ", object@dtime, "\n",
                  sep="")
              cat("  Number of Samples    : ", length(object@time), "\n",
                  sep="")
              cat("  Sampling Begins      : ",
                  paste(object@time[1]), "\n", sep="")
              cat("  Sampling Ends        : ",
                  paste(object@time[length(object@time)]), "\n", sep="")
              cat("  Total Duration (d)   : ",
                  difftime(trange[2], trange[1], units="days"), "\n", sep="")
              drange <- range(object@depth, na.rm=TRUE)
              cat("  Measured depth range : [",
                  drange[1], ", ", drange[2], "]\n", sep="")
              if (length(names(object@concurrentData)) > 0) {
                  cat("  Other variables      : ",
                      names(object@concurrentData), "\n")
              }
          })

setMethod("plotTDR", signature(x="POSIXt", y="numeric"),
          function(x, y, concurVars=NULL, xlim=NULL, depth.lim=NULL,
                   xlab="time (dd-mmm hh:mm)", ylab.depth="depth (m)",
                   concurVarTitles=deparse(substitute(concurVars)),
                   xlab.format="%d-%b %H:%M", sunrise.time="06:00:00",
                   sunset.time="18:00:00", night.col="gray60", dry.time=NULL,
                   phase.factor=NULL, plot.points=FALSE, interact=TRUE,
                   key=TRUE, cex.pts=0.4, ...) {
              stopifnot(identical(length(x), length(y)), is.vector(y))
              .plotTDR(time=x, depth=y, concurVars=concurVars,
                       xlim=xlim, depth.lim=depth.lim, xlab=xlab,
                       ylab.depth=ylab.depth,
                       concurVarTitles=concurVarTitles,
                       xlab.format=xlab.format,
                       sunrise.time=sunrise.time,
                       sunset.time=sunset.time,
                       night.col=night.col, dry.time=dry.time,
                       phase.factor=phase.factor,
                       interact=interact, key=key,
                       cex.pts=cex.pts, ...)
          })

setMethod("plotTDR", signature(x="TDR", y="missing"),
          function(x, y, concurVars, concurVarTitles, ...) {
              if (!missing(concurVars)) {
                  ccd <- getCCData(x, concurVars)
                  if (!missing(concurVarTitles)) {
                      lcvt <- length(concurVarTitles)
                      lcv <- length(concurVars)
                      stopifnot(identical(lcvt, lcv))
                  } else concurVarTitles <- colnames(ccd)
              } else if (missing(concurVars) && missing(concurVarTitles)) {
                  ccd <- concurVarTitles <- NULL
              }
              .plotTDR(time=getTime(x), depth=getDepth(x),
                       concurVars=ccd,
                       concurVarTitles=concurVarTitles, ...)
          })

###_  . TDRcalibrate
setMethod("show", signature=signature(object="TDRcalibrate"),
          definition=function(object) {
              mCall <- gsub(" = ", "=", gsub("^ +", "", deparse(object@call)))
              dry <- object@gross.activity$activity == "L"
              dd <- length(unique(object@gross.activity$ phase.id[dry]))
              wet <- object@gross.activity$activity == "W"
              wetz <- object@gross.activity$activity == "Z"
              ww <- length(unique(object@gross.activity$ phase.id[wet | wetz]))
              cat("Depth calibration -- Class", class(object), "object\n")
              cat("  Call                              : ", mCall, "\n", sep="")
              cat("  Source file                       : ", object@tdr@file, "\n",
                  sep="")
              cat("  Containing TDR of class           : ", class(object@tdr),
                  "\n", sep="")
              cat("  Number of dry phases              : ", dd, "\n", sep="")
              cat("  Number of aquatic phases          : ", ww, "\n", sep="")
              cat("  Number of dives detected          : ",
                  max(object@dive.activity$dive.id, na.rm=TRUE), "\n", sep="")
              cat("  Dry threshold used (s)            : ", object@dry.thr, "\n",
                  sep="")
              cat("  Aquatic theshold used (s)         : ", object@wet.thr, "\n",
                  sep="")
              cat("  Dive threshold used (depth units) : ", object@dive.thr,
                  sep="")
              if (is(object@tdr, "TDRspeed")) {
                  cat("\n  Speed calibration coefficients    : a=",
                      format(object@speed.calib.coefs[1], digits=2), "; b=",
                      format(object@speed.calib.coefs[2], digits=2), "\n",
                      sep="")
              } else cat("\n", sep="")
          })

".plotTDRcalibratePhases" <- function(x, diveNo=seq(max(getDAct(x, "dive.id"))),
                                      concurVars, surface=FALSE, ...)
{
    if (!is(x, "TDRcalibrate")) stop("x must be a TDRcalibrate object")
    ell <- list(...)
    diveNo <- sort(diveNo)
    diveids <- getDAct(x, "dive.id")
    tdr <- getTDR(x)
    if (max(unique(diveids)) < 1) {
        ok <- seq(along=slot(tdr, "depth"))
    } else if (surface) {
        dives <- diveids %in% diveNo
        postdiveids <- getDAct(x, "postdive.id")
        postdives <- postdiveids %in% diveNo
        ok <- which(dives | postdives)
    } else ok <- .diveIndices(diveids, diveNo)
    newtdr <- tdr[ok]
    alltimes <- getTime(tdr)
    newtimes <- getTime(newtdr)
    times.ok <- alltimes >= newtimes[1] & alltimes <= newtimes[length(newtimes)]
    fulltimes <- alltimes[times.ok]
    labs <- getDPhaseLab(x)[ok]
    drys <- getGAct(x, "activity")[times.ok]
    drys[drys == "Z"] <- "L"; drys <- drys[, drop=TRUE]
    dry.time <- fulltimes[drys == "L"]
    ell$x <- newtdr
    ell$phase.factor <- labs
    if (length(dry.time) > 0L) ell$dry.time <- dry.time
    if (!missing(concurVars)) {
        if (!is.character(concurVars))
            stop("concurVars must be of class character")
        ell$concurVars <- concurVars
    }
    do.call(plotTDR, args=ell)
}
setMethod("plotTDR", signature(x="TDRcalibrate", y="missing"),
          function(x, y, what=c("phases", "dive.model"),
                   diveNo=seq(max(getDAct(x, "dive.id"))), ...) {
              what <- match.arg(what)
              switch(what,
                     phases = {
                         .plotTDRcalibratePhases(x, diveNo=diveNo, ...)
                     },
                     dive.model = { plotDiveModel(x, diveNo=diveNo) })
          })

###_  . diveModel
setMethod("show", signature=signature(object="diveModel"),
          definition=function(object) {
              ## Lots stolen from print.smooth.spline()
              digits <- getOption("digits")
              cat("Dive model -- Class",
                  class(object), "object\n")
              if(!is.null(cl <- object@dive.spline$call)) {
                  cat("Call:\n")
                  dput(cl, control=NULL)
              }
              ip <- object@dive.spline$iparms
              cv <- cl$cv
              if(is.null(cv)) cv <- FALSE else if(is.name(cv)) cv <- eval(cv)
              cat("\nSmoothing Parameter  spar=",
                  format(object@dive.spline$spar, digits=digits),
                  " lambda=", format(object@dive.spline$lambda, digits=digits),
                  if (ip["ispar"] != 1L) {
                      paste("(", ip["iter"], " iterations)", sep="")
                      }, "\n", sep="")
              cat("Equivalent Degrees of Freedom :",
                  format(object@dive.spline$df, digits=digits), "\n")
              cat("Penalized Criterion           :",
                  format(object@dive.spline$pen.crit, digits=digits), "\n")
              cat(ifelse(cv,
                         "PRESS                         : ",
                         "GCV                           : "),
                  format(object@dive.spline$cv.crit, digits=digits),
                  "\n", sep="")
              cat("Observed N                    : ",
                  nrow(object@label.matrix), "\n", sep="")
              cat("Modelled N                    : ",
                  length(object@dive.spline$data$x), "\n", sep="")
              cat("Modelled N (distinct)         : ",
                  length(object@dive.spline$x), "\n", sep="")
              cat("Derivative evaluated at       : ",
                  length(object@spline.deriv$x), " points", "\n", sep="")
              cat("Descent ends after            : ",
                  object@descent.crit, " steps in model", "\n", sep="")
              cat("Ascent begins after           : ",
                  object@ascent.crit, " steps in model", "\n", sep="")
              cat("Descent critical rate         : ",
                  object@descent.crit.rate, "\n", sep="")
              cat("Ascent critical rate          : ",
                  object@ascent.crit.rate, "\n", sep="")
          })

setMethod("plotDiveModel", signature(x="diveModel", y="missing"),
          function(x, diveNo) {
              if (missing(diveNo)) diveNo <- "Unknown"
              diveM <- x
              times <- diveM@dive.spline$data$x
              depths <- diveM@dive.spline$data$y
              times.s <- diveM@dive.spline$x
              depths.s <- diveM@dive.spline$y
              times.deriv <- diveM@spline.deriv$x
              depths.deriv <- diveM@spline.deriv$y
              d.crit <- diveM@descent.crit
              a.crit <- diveM@ascent.crit
              d.crit.rate <- diveM@descent.crit.rate
              a.crit.rate <- diveM@ascent.crit.rate
              plotDiveModel(x=times, y=depths, times.s=times.s,
                            depths.s=depths.s, d.crit=d.crit, a.crit=a.crit,
                            times.deriv=times.deriv,
                            depths.deriv=depths.deriv, diveNo=diveNo,
                            d.crit.rate=d.crit.rate, a.crit.rate=a.crit.rate)
          })

setMethod("plotDiveModel",
          signature(x="TDRcalibrate", y="missing"),
          function(x, diveNo) {
              if (length(diveNo) != 1L)
                  stop("Only one dive's phase model can be plotted")
              diveM <- getDiveModel(x, diveNo)
              dive <- extractDive(x, diveNo)
              times <- getTime(dive)
              times <- as.numeric(times - times[1])
              depths <- getDepth(dive)
              times.s <- diveM@dive.spline$x
              depths.s <- diveM@dive.spline$y
              times.deriv <- diveM@spline.deriv$x
              if (length(times) < 4L) {
                  ff <- times[length(times)] / times.s[length(times.s)]
                  times.s <- times.s * ff
                  times.deriv <- times.deriv * ff
              }
              depths.deriv <- diveM@spline.deriv$y
              d.crit <- diveM@descent.crit
              a.crit <- diveM@ascent.crit
              d.crit.rate <- diveM@descent.crit.rate
              a.crit.rate <- diveM@ascent.crit.rate
              plotDiveModel(x=times, y=depths, times.s=times.s,
                            depths.s=depths.s, d.crit=d.crit, a.crit=a.crit,
                            diveNo=diveNo, times.deriv=times.deriv,
                            depths.deriv=depths.deriv,
                            d.crit.rate=d.crit.rate, a.crit.rate=a.crit.rate)
          })

setMethod("plotDiveModel",
          signature(x="numeric", y="numeric"),
          function(x, y, times.s, depths.s, d.crit, a.crit, diveNo=1,
                   times.deriv, depths.deriv, d.crit.rate, a.crit.rate) {
              times <- x
              depths <- -abs(y)
              depths.s <- -abs(depths.s)
              descent.c1 <- times.deriv < times[d.crit]
              descent.c2 <- depths.deriv > d.crit.rate
              descent <- descent.c1 & descent.c2
              ascent.c1 <- times.deriv > times[a.crit]
              ascent.c2 <- depths.deriv < a.crit.rate
              ascent <- ascent.c1 & ascent.c2
              layout(matrix(1:2, ncol=1))
              old.par <- par(no.readonly=TRUE)
              on.exit(par(old.par))
              par(mar=c(3, 4, 0, 1) + 0.1, las=1)
              plot(times, depths, type="o", axes=FALSE, pch=19, cex=0.5,
                   frame.plot=TRUE, ylab="Depth",
                   ylim=range(depths, depths.s, na.rm=TRUE))
              axis(side=1)
              axis(side=2, at=pretty(c(depths, depths.s)),
                   labels=rev(pretty(-c(depths, depths.s))), las=1)
              lines(times.s, depths.s, lty=2, col="green")
              lines(times[seq(d.crit)], depths[seq(d.crit)], col="blue")
              lines(times[a.crit:length(x)],
                    depths[a.crit:length(x)], col="lightblue")
              legend("top", ncol=2, title=paste("Dive:", diveNo),
                     legend=c("original", "smoothed",
                       "descent", "ascent"), lty=c(1, 2, 1, 1),
                     col=c("black", "green",
                       "blue", "lightblue"), cex=0.7)
              plot(times.deriv, depths.deriv, xlab="Time index",
                   ylab="First derivative", type="l", cex=0.3)
              points(times.deriv[descent], depths.deriv[descent],
                     col="blue", cex=0.3)
              points(times.deriv[ascent], depths.deriv[ascent],
                     col="lightblue", cex=0.3)
              abline(h=c(d.crit.rate, a.crit.rate),
                     v=c(times[d.crit], times[a.crit]), lty=2)
              text(2, c(d.crit.rate, a.crit.rate),
                   labels=c(expression(paste("descent ", hat(q))),
                     expression(paste("ascent ", hat(q)))),
                   pos=c(3, 1), cex=0.7)
              text(c(times[d.crit], times[a.crit]), 0,
                   labels=c("descent", "ascent"), pos=1, cex=0.7)
          })

###_  . plotBouts
setMethod("plotBouts", signature(fit="nls"),
          function(fit, ...) {
              ncoefs <- as.character(length(coef(fit)))
              if (! (ncoefs == "4" || ncoefs == "6")) {
                  msg <- paste("fitted model must have 4 (2-process) or",
                               "6 (3-process) coefficients")
                  stop(msg)
              }
              switch(ncoefs,
                     "4" = {
                         plotBouts2.nls(fit=fit,
                                        lnfreq=eval.parent(fit$data), ...)
                     },
                     "6" = {
                         plotBouts3.nls(fit=fit,
                                        lnfreq=eval.parent(fit$data), ...)
                     })
          })
setMethod("plotBouts", signature(fit="mle"),
          function(fit, x, ...) {
              ncoefs <- as.character(length(coef(fit)))
              if (! (ncoefs == "3" || ncoefs == "5")) {
                  msg <- paste("fitted model must have 3 (2-process) or",
                               "5 (3-process) coefficients")
                  stop(msg)
              }
              switch(ncoefs,
                     "3" = {
                         plotBouts2.mle(fit=fit, x=x, ...)
                     },
                     "5" = {
                         stop("To be implemented")
                     })
          })

###_  . plotZOC
setMethod("plotZOC", signature(x="TDR", y="matrix"),
          function(x, y, xlim, ylim, ylab="Depth (m)", ...) {
              .plotZOCfilters(x=x, zoc.filter=y, xlim=xlim, ylim=ylim,
                              ylab=ylab, ...)
          })

setMethod("plotZOC", signature(x="TDR", y="TDRcalibrate"),
          function(x, y, xlim, ylim, ylab="Depth (m)", ...) {
              .plotZOCtdrs(x=x, y=y, xlim=xlim, ylim=ylim, ylab=ylab, ...)
          })


###_ + Accessors

###_  . TDR and TDRspeed
setMethod("getFileName", signature(x="TDR"), function(x) x@file)

if (getRversion() < "2.11.0") {
    .POSIXct <- function(xx, tz=NULL) {
        structure(xx, class=c("POSIXt", "POSIXct"), tzone=tz)
    }
}
setMethod("getTime", signature(x="TDR"),
          function(x) {
              xx <- x@time
              if (getRversion() >= "2.12.0") {
                  .POSIXct(unclass(xx), attr(xx, "tzone"))
              } else xx
          })

setMethod("getDepth", signature(x="TDR"), function(x) x@depth)

".speedCol" <- function(x)
{
    ## Value: column number where speed is located in x
    ## --------------------------------------------------------------------
    ## Arguments: x=data frame
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    dataNames <- names(x)
    colN <- dataNames %in% .speedNames
    if (length(which(colN)) != 1)
        stop("the column number for speed could not be determined")
    which(colN)
}
setMethod("getSpeed", signature(x="TDRspeed"), function(x) {
    ccData <- x@concurrentData
    speedCol <- .speedCol(ccData)
    ccData[, speedCol]
})

setMethod("getDtime", signature(x="TDR"), function(x) x@dtime)

## Get entire data frame
setMethod("getCCData", signature(x="TDR", y="missing"), function(x) {
    if (nrow(x@concurrentData) > 0) {
        x@concurrentData
    } else stop("No concurrent data are available")
})
## Get named component(s) of the data frame
setMethod("getCCData", signature(x="TDR", y="character"), function(x, y) {
    if (nrow(x@concurrentData) < 1) stop("No concurrent data are available")
    ccd <- getCCData(x)
    ccdnames <- names(ccd)
    ok <- ccdnames %in% y
    bady <- !y %in% ccdnames
    if (length(which(ok)) < 1) {
        stop(paste("y must contain at least one of the names of",
                   "the concurrent data frame"))
    } else if (any(bady)) {
        warning("components: ", y[bady], " could not be found and were ignored")
    }
    ccdf <- as.data.frame(ccd[, ok])
    names(ccdf) <- ccdnames[ok]
    ccdf
})

###_  . TDRcalibrate
setMethod("getTDR", signature(x="TDRcalibrate"), function(x) x@tdr)

## access the entire list
setMethod("getGAct", signature(x="TDRcalibrate", y="missing"),
          function(x) x@gross.activity)
## access only a named element
setMethod("getGAct", signature(x="TDRcalibrate", y="character"),
          function(x, y) x@gross.activity[[y]])

## access entire data frame
setMethod("getDAct", signature(x="TDRcalibrate", y="missing"),
          function(x) x@dive.activity)
## access only a certain column
setMethod("getDAct", signature(x="TDRcalibrate", y="character"),
          function(x, y) x@dive.activity[, y])

## access the entire factor
setMethod("getDPhaseLab", signature(x="TDRcalibrate", diveNo="missing"),
          function(x) x@dive.phases)
## access only those from certain dives
setMethod("getDPhaseLab", signature(x="TDRcalibrate", diveNo="numeric"),
          function(x, diveNo) {
              ctdr <- getTDR(x)
              phases <- x@dive.phases
              okpts <- .diveIndices(getDAct(x, "dive.id"), diveNo)
              phases[okpts]
          })

## access the entire object
setMethod("getDiveModel", signature(x="TDRcalibrate", diveNo="missing"),
          function(x) x@dive.models)
## access only those from certain dives -- simplify if only one
setMethod("getDiveModel", signature(x="TDRcalibrate", diveNo="numeric"),
          function(x, diveNo) {
              dml <- x@dive.models
              tryCatch({
                  ok <- .diveMatches(names(dml), diveNo)
                  diveNo.ok <- diveNo[ok]
                  dm <- x@dive.models[diveNo.ok]
                  if (length(diveNo.ok) == 1L) dm[[1]] else dm
              })
          })

## Basic diveModel
setMethod("getDiveDeriv", signature(x="diveModel"),
          function(x, phase=c("all", "descent", "bottom", "ascent")) {
              phase <- match.arg(phase)
              d.crit <- x@descent.crit
              a.crit <- x@ascent.crit
              switch(phase,
                     all = {x@spline.deriv},
                     descent = {
                         spd <- x@spline.deriv
                         t.crit <- x@dive.spline$data$x[d.crit]
                         descent <- which(spd$x < t.crit)
                         spd$x <- spd$x[descent]
                         spd$y <- spd$y[descent]
                         spd
                     },
                     bottom = {
                         spd <- x@spline.deriv
                         t.desc.crit <- x@dive.spline$data$x[d.crit]
                         t.asc.crit <- x@dive.spline$data$x[a.crit]
                         bottom <- which(spd$x >= t.desc.crit &
                                         spd$x <= t.asc.crit)
                         spd$x <- spd$x[bottom]
                         spd$y <- spd$y[bottom]
                         spd
                     },
                     ascent = {
                         spd <- x@spline.deriv
                         t.crit <- x@dive.spline$data$x[a.crit]
                         ascent <- which(spd$x > t.crit)
                         spd$x <- spd$x[ascent]
                         spd$y <- spd$y[ascent]
                         spd
                     })
          })
## TDRcalibrate -- do all dives or selection.  Simplify if only one
setMethod("getDiveDeriv", signature(x="TDRcalibrate"),
          function(x, diveNo, phase=c("all", "descent", "bottom", "ascent")) {
              if (missing(diveNo)) diveNo <- seq(max(getDAct(x, "dive.id")))
              phase <- match.arg(phase)
              dl <- lapply(diveNo, function(k) {
                  dm <- getDiveModel(x, diveNo=k)
                  getDiveDeriv(dm, phase=phase)
              })
              names(dl) <- diveNo
              if (length(diveNo) == 1L) dl[[1]] else dl
          })

setMethod("getSpeedCoef", signature(x="TDRcalibrate"),
          function(x) x@speed.calib.coefs)


###_ + Coercions and Replacements
setAs("TDR", "data.frame", function(from) {
    file.src <- from@file
    dtime <- from@dtime
    val <- data.frame(time=from@time, depth=from@depth, getCCData(from))
    attr(val, "file") <- file.src
    attr(val, "dtime") <- dtime
    val
})
setMethod("as.data.frame", signature("TDR"),
          function(x, row.names=NULL, optional=FALSE) {
              as(x, "data.frame")
          })

setAs("TDR", "TDRspeed", function(from) {
    new("TDRspeed", file=from@file, time=from@time, depth=from@depth,
        dtime=from@dtime, concurrentData=from@concurrentData)
})
setMethod("as.TDRspeed", signature("TDR"), function(x) as(x, "TDRspeed"))

setReplaceMethod("depth", signature(x="TDR", value="numeric"),
                 function(x, value) {
                     orig <- getDepth(x)
                     if (length(orig) != length(value))
                         stop(paste("replacement must have length:",
                                    length(orig),
                                    "(i.e. same length as original)"))
                     x@depth <- value
                     x
                 })

setReplaceMethod("speed", signature(x="TDRspeed", value="numeric"),
                 function(x, value) {
                     ccData <- x@concurrentData
                     speedCol <- .speedCol(ccData)
                     if (length(ccData[, speedCol]) != length(value))
                         stop(paste("replacement must have length:",
                                    length(ccData[, speedCol]),
                                    "(i.e. same length as original)"))
                     x@concurrentData[, speedCol] <- value
                     x
                 })

setReplaceMethod("ccData", signature(x="TDR", value="data.frame"),
                 function(x, value) {
                     ccDataN <- nrow(getCCData(x))
                     valueN <- nrow(value)
                     if (ccDataN != valueN)
                         stop(paste("replacement must have:", ccDataN,
                                    "rows (i.e. same as original)"))
                     x@concurrentData <- value
                     x
                 })

###_ + Subsetting
setMethod("[", signature(x="TDR", i="numeric", j="missing", drop="missing"),
          function(x, i, j, ..., drop) {
    new(class(x), file=getFileName(x), dtime=getDtime(x), time=getTime(x)[i],
        depth=getDepth(x)[i],
        concurrentData=tryCatch(getCCData(x)[i, , drop=FALSE],
          error=function(k) data.frame()))
})


###_ + Generators and Summaries
"createTDR" <- function(time, depth,
                        concurrentData=data.frame(matrix(ncol=0,
                          nrow=length(time))),
                        speed=FALSE, dtime, file)
{
    ## Value: An object of TDR or TDRspeed class.  Useful to recreate
    ## objects once depth has been zoc'ed and speed calibrated for further
    ## analyses.
    ## --------------------------------------------------------------------
    ## Arguments: see class definitions
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    if (missing(dtime)) dtime <- .getInterval(time)
    if(speed) {
        new("TDRspeed", time=time, depth=depth, concurrentData=concurrentData,
            dtime=dtime, file=file)
    } else {
        new("TDR", time=time, depth=depth, concurrentData=concurrentData,
            dtime=dtime, file=file)
    }
}

setMethod("extractDive", signature(obj="TDR", diveNo="numeric",
                                   id="numeric"), # for TDR object
          function(obj, diveNo, id) {
              if (length(id) != length(getTime(obj))) {
                  stop ("id and obj must have equal number of rows")
              }
              okpts <- .diveIndices(id, unique(diveNo))
              if (is(obj, "TDRspeed")) {
                  new("TDRspeed", time=getTime(obj)[okpts],
                      depth=getDepth(obj)[okpts],
                      concurrentData=getCCData(obj)[okpts, , drop=FALSE],
                      dtime=getDtime(obj), file=obj@file)
              } else {
                  new("TDR", time=getTime(obj)[okpts],
                      depth=getDepth(obj)[okpts],
                      concurrentData=getCCData(obj)[okpts, , drop=FALSE],
                      dtime=getDtime(obj), file=obj@file)
              }
          })

setMethod("extractDive",                # for TDRcalibrate
          signature(obj="TDRcalibrate", diveNo="numeric", id="missing"),
          function(obj, diveNo) {
              ctdr <- getTDR(obj)
              okpts <- .diveIndices(getDAct(obj, "dive.id"),
                                    unique(diveNo))
              if (is(ctdr, "TDRspeed")) {
                  new("TDRspeed", time=getTime(ctdr)[okpts],
                      depth=getDepth(ctdr)[okpts],
                      concurrentData=getCCData(ctdr)[okpts, , drop=FALSE],
                      dtime=getDtime(ctdr), file=ctdr@file)
              } else {
                  new("TDR", time=getTime(ctdr)[okpts],
                      depth=getDepth(ctdr)[okpts],
                      concurrentData=getCCData(ctdr)[okpts, , drop=FALSE],
                      dtime=getDtime(ctdr), file=ctdr@file)
              }
          })

setMethod("timeBudget",            # a table of general attendance pattern
          signature(obj="TDRcalibrate", ignoreZ="logical"),
          function(obj, ignoreZ) {
              act <- getGAct(obj, "activity")
              tt <- getTime(getTDR(obj))
              interval <- getDtime(getTDR(obj))
              if (ignoreZ) {            # ignore the short baths
                  act[act == "Z"] <- "L"
                  attlist <- .rleActivity(tt, act, interval)
                  actlabel <- rle(as.vector(act))$values
                  phase.no <- seq(along=actlabel)
              } else {                  # count the short baths
                  attlist <- getGAct(obj)
                  actlabel <- rle(as.vector(act))$values
                  phase.no <- seq(along=actlabel)
              }
              data.frame(phase.no=phase.no, activity=actlabel,
                         beg=attlist[[3]], end=attlist[[4]],
                         row.names=NULL)
          })


###_ + Methods for bec2 and bec3 are in bouts.R

## This is to avoid Collate issues in DESCRIPTION


###_ + Emacs local variables
## Local variables:
## allout-layout: (+ : 0)
## End:
