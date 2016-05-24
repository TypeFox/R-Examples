#
# 	dive.R
#
#	$Revision: 1.30 $	$Date: 2014/01/28 02:16:00 $
#
###################################################################
#  
# Dive profiles
#
# e.g.
#      dive(c(21,30),c(5,5),c(0,120),c(18,30),c(5,5))
#	  means: 21 m for 30 min, 5@5, 2 hour surface, 18 m for 30 min, 5@5
#
#      dive(nitrox(0.32), c(25,30), 5, nitrox(0.5), c(5,5))
#                 means dive on EANx 32 to 25 metres for 30 min,
#                 ascending to 5 metres, switching to EANx 50,
#                 5@5 then surfacing.
#

dive <- function(..., begin=0, end=0, tanklist=NULL) {
  X <- list(...)
  if(length(X) == 0)
    stop("No dive profile supplied")

  if(!is.na(begin))
    stopifnot(is.numeric(begin) && length(begin) == 1)
  if(!is.na(end))
    stopifnot(is.numeric(end) && length(end) == 1)

  tanks.given <- !is.null(tanklist)
  gases.given <- tanks.given || any(unlist(lapply(X, is.gas)))
  tanks.implied <- FALSE

  # initialise current state and rules
  state <- data.frame(time=0, depth=begin,
                      fO2=air$fO2, fN2=air$fN2, fHe=air$fHe, tankid=1)
  rate.down <- descent(30)  # default descent rate
  rate.up <- ascent(18)    # default ascent rate
  if(!tanks.given)
    tanklist <- list(air)  # default tank
  
  # set up data frame to record the entire dive profile
  # NB: gas fractions apply to the time interval [i,i+1]
  df <- state
  now <- 1

  # Some entries of X may be labelled (tag=value)
  tagged <- nzchar(names(X))
  if(length(tagged) == 0) tagged <- rep(FALSE, length(X)) 
  
  change <- function(state, duration=NULL, depth=NULL,
                     fO2=NULL, fN2=NULL, fHe=NULL, tankid=NULL) {
    if(!is.null(duration)) state$time <- state$time + duration
    if(!is.null(depth))    state$depth <- depth
    if(!is.null(fO2))      state$fO2 <- fO2
    if(!is.null(fN2))      state$fN2 <- fN2
    if(!is.null(fHe))      state$fHe <- fHe
    if(!is.null(tankid))   state$tankid <- tankid
    return(state)
  }

  gasequal <- function(g1, g2) { identical(all.equal(g1, g2), TRUE) }

  # examine arguments one-by-one
  for(i in seq(X)) {
    Y <- X[[i]]
    # ....................................................................
    if(is.rate(Y)) # change the ascent or descent rate
      if(Y$up) rate.up <- Y else rate.down <- Y
    else if(is.gas(Y)) {
      # switch to new gas for next time interval
      if(tanks.given)
        stop(paste("If tanklist is given, gas switches should indicate",
                   "which tank is used (tank=number or tank=name)"))
      matched <- unlist(lapply(tanklist, gasequal, g2=Y))
      if(any(matched)) 
        tankid <- min(which(matched))
      else {
        # new gas
        if(now == 1) 
          tanklist <- list(Y)
        else 
          tanklist <- append(tanklist, list(Y))
        tankid <- length(tanklist)
      }
      state <- change(state, fO2=Y$fO2, fN2=Y$fN2, fHe=Y$fHe, tankid=tankid)
      df[now, ] <- state
    } else if(tanks.given && tagged[i] && names(X)[i] == "tank") {
      # argument tank=name or tank=number
      # switch to new tank Y for next time interval
      if(is.null(tanklist))
        stop("Cannot change tank: no tanks specified")
      Z <- tanklist[[Y]]
      state <- change(state, fO2=Z$fO2, fN2=Z$fN2, fHe=Z$fHe, tankid=Y)
      df[now, ] <- state
    } else if(is.vector(Y)) {
      # waypoint
      if(length(Y) > 2)
        stop("Vector of length > 2 is not a recognised format")
      newdepth <- Y[1]
      currentdepth <- state$depth
      if(is.na(currentdepth)) {
        # re-initialise depth
        state <- change(state, depth=newdepth)
      } else if(newdepth != currentdepth) {
        # Ascend or descend to this new depth (making new waypoint)
        duration <- timetaken(currentdepth, newdepth, rate.up, rate.down)
        state <- change(state, duration, newdepth)
        df <- rbind(df, state)
        now <- now + 1
      }
      # Now stay at this depth for the specified time (making new waypoint)
      if(length(Y) > 1) {
        duration <- Y[2]
        state <- change(state, duration)
        df <- rbind(df, state)
        now <- now + 1
      }
    } else if(is.data.frame(Y)) {
      # dive profile data
      if(ncol(Y) != 2)
        stop("Data frame should have 2 columns")
      newtimes <- Y[,1]
      newdepths <- Y[,2]
      nsteps <- nrow(Y)

      if(inherits(newtimes, "difftime")) {
        # convert to minutes
        units(newtimes) <- "mins"
        newtimes <- as.numeric(newtimes)
      } else if(is.numeric(newtimes)) {
        # assume seconds -- convert to minutes
        message("Elapsed times are assumed to be given in seconds (and are converted to minutes)")
        newtimes <- newtimes/60
      } else if(is.character(newtimes)) {
        # assume mm:ss format
        zip <- function(x){
          x <- as.numeric(x)
          n <- length(x)
          secs <- sum(x * 60^((n-1):0))
          return(secs/60)
        }
        newtimes <- unlist(lapply(strsplit(newtimes, ":"), zip))
      } else 
        stop("Unknown format for elapsed times")

      if(any(is.na(newtimes)))
        stop("Unable to convert times")
      if(!all(diff(newtimes) > 0))
        stop("Elapsed times are not an increasing sequence")

      # convert elapsed times
      newtimes <- newtimes + state$time
      
      # make data frame 
      newdf <- data.frame(time=newtimes,
                          depth=newdepths,
                          fO2=state$fO2,
                          fN2=state$fN2,
                          fHe=state$fHe)
      if(!is.null(tanklist))
        newdf <- cbind(newdf, tankid=state$tankid)

      currentdepth <- state$depth
      if(is.na(currentdepth)) {
        # dive not yet started: take X as dive data
        df <- newdf
      } else {
        # dive already started - transition to X
        if(newdepths[1] != currentdepth) {
          # First ascend/descend to starting depth (on current gas)
          duration <- timetaken(currentdepth, newdepths[1], rate.up, rate.down)
          newdf$time <- newdf$time + duration
        } else if(newtimes[1] == state$time) {
          # last row of df duplicates first row of newdf
          df <- df[-now, ]
        }         
        # tack on the data frame
        df <- rbind(df, newdf)
      }
      now <- nrow(df)
      state <- df[now, ]
      
    } else if(is.dive(Y)) {
      #### dive object
      Ydata <- Y$data
      # First ascend/descend to depth (on current gas)
      newdepth <- Ydata$depth[1]
      if(newdepth != state$depth) {
        duration <- timetaken(state$depth, newdepth, rate.up, rate.down)
        state <- change(state, duration, newdepth)
        df <- rbind(df, state)
        now <- now + 1
      } else {
          # last row of df duplicates first row of newdf
          df <- df[-now, ]
      }         
      # Now determine breathing gases and alter 'Ydata'
      if(gases.given) {
        # The arguments in the call to dive() specify some gases.
        # Override the gases used in Y.
        # Assume Y is conducted using current gas
        Ydata$tankid <- state$tankid
        Ydata$fO2    <- state$fO2
        Ydata$fN2    <- state$fN2
        Ydata$fHe    <- state$fHe
      } else {
        # No gases were specified in the call to dive().
        if(now == 1) {
          # We haven't started diving yet, so no tanks have been breathed
          # and the tank list is so far undetermined.
          # Copy the (explicit or implicit) tank list from Y
          # Copy tank selection from Y
          tanklist <- Y$tanklist
          tanks.implied <- tanks.implied || Y$tanks.given
        } else if(!Y$tanks.given) {
          # Y does not have any user-defined tank information;
          # override gases in Y
          # Assume Y is conducted using current gas
          Ydata$tankid <- state$tankid
          Ydata$fO2    <- state$fO2
          Ydata$fN2    <- state$fN2
          Ydata$fHe    <- state$fHe
        } else {
          # Y has user-defined tank information;
          # our dive has already started on some tank;
          # merge tank lists
          tanks.implied <- TRUE
          tid <- df$tankid
          Ytid <- Ydata$tankid
          if(is.factor(tid) && is.factor(Ytid)) {
            # tanks identified by factor/character labels
            if(!identical(tanklist, Y$tanklist)
               || !identical(levels(tid), levels(Ytid))) {
              # Factors are not compatible
              # Take union of factor levels
              tidN <- length(levels(tid))
              YtidN <- length(levels(Ytid))
              lev <- c(levels(tid), levels(Ytid))
              tid <- factor(tid, levels=lev)
              Ytid <- factor(as.integer(Ytid) + tidN, levels=seq(tidN + YtidN))
              levels(Ytid) <- lev
              df$tankid <- tid
              Ydata$tankid <- Ytid
              tanklist <- append(tanklist, Y$tanklist)
            }
          } else {
            # tanks will be numbered
            df$tankid <- as.integer(tid)
            Ydata$tankid <- as.integer(Ydata$tankid) + length(tanklist)
            tanklist <- append(tanklist, Y$tanklist)
          }
          # update gas fractions in 'Ydata'
          Ydata <- reconcile.df(Ydata, tanklist)
        }
      }
      # Then tack on the new dive
      Ydata$time <- (Ydata$time - Ydata$time[1]) + state$time
      df <- rbind(df, Ydata)
      state <- df[nrow(df), ]
      now <- nrow(df) + 1
    } else 
    stop(paste("The format of argument", i, "is not recognised\n"))
    # ....................................................................
  }
  # End of loop - all arguments examined.
  #
  if(!is.na(end) && state$depth != end) {
    # Don't forget to surface.. (or whatever)
    duration <- timetaken(state$depth, end, rate.up, rate.down)
    state <- change(state, duration, end)
    df <- rbind(df, state)
    now <- now + 1
  }
  # check dive is all under water!
  if(any(df$depth < 0)) {
    warning("Some depths were less than 0; they were reset to 0")
    df$depth <- pmax(0, df$depth)
  }
  # was it an air dive?
  airdive <- all(with(df, gasnames(fO2, fN2, fHe)) == "air")

  # convert character tank names to factor levels
  # so that they can always be coerced to integer
  if(is.character(df$tankid)) {
    nama <- names(tanklist)
    if(is.null(nama))
      nama <- paste(seq(tanklist))
    df$tankid <- factor(df$tankid, levels=nama)
  }
  
  # dive depth/time profile
  if(nrow(df) <= 1)
    stop("Need at least two waypoints for a dive")
  pro <- approxfun(df$time, df$depth, method="linear", rule=2, ties="ordered")
  # 
  result <- list(profile=pro, data=df, airdive=airdive, tanklist=tanklist,
                 tanks.given=tanks.given || tanks.implied)
  class(result) <- c("dive", class(result))
  return(result)
}

is.dive <- function(x) { inherits(x,"dive")}

####################################################################
# Entries in dive objects

depths.dive <- function(d) {
  stopifnot(is.dive(d))
  d$data$depth
}

times.dive <- function(d) {
  stopifnot(is.dive(d))
  d$data$time
}

durations.dive <- function(d) {
  stopifnot(is.dive(d))
  diff(times.dive(d))
}

"depths.dive<-" <- function(d, value) {
  stopifnot(is.dive(d))
  stopifnot(is.vector(value) && is.numeric(value))
  if(length(value) != nrow(d$data))
    stop("incorrect length of replacement vector for depths")
  d$data$depth <- value
  # update dive depth/time profile
  d$profilepro <- approxfun(d$data$time, d$data$depth,
                            method="linear", rule=2, ties="ordered")
  # 
  return(d)
}

"times.dive<-" <- function(d, value) {
  stopifnot(is.dive(d))
  stopifnot(is.vector(value) && is.numeric(value))
  if(length(value) != nrow(d$data))
    stop("incorrect length of replacement vector for times")
  d$data$time <- value
  # update dive depth/time profile
  d$profilepro <- approxfun(d$data$time, d$data$depth,
                            method="linear", rule=2, ties="ordered")
  return(d)
}

"durations.dive<-" <- function(d, value) {
  stopifnot(is.dive(d))
  stopifnot(is.vector(value) && is.numeric(value))
  if(length(value) != nrow(d$data) - 1)
    stop("incorrect length of replacement vector for durations")
  d$data$time <- d$data$time[1] + cumsum(c(0, value))
  # update dive depth/time profile
  d$profilepro <- approxfun(d$data$time, d$data$depth,
                            method="linear", rule=2, ties="ordered")
  return(d)
}

whichtank <- function(d) {
  stopifnot(is.dive(d))
  return(d$data$tankid)
}

"whichtank<-" <- function(d, value) {
  stopifnot(is.dive(d))
  DF <- d$data
  TL <- d$tanklist
  ndf <- nrow(DF)

  checkcollapse <- TRUE
  if(length(value) == 1) { 
    value <- rep(value, ndf)
    checkcollapse <- FALSE
  } else if(length(value) != ndf) 
    stop("Replacement value has the wrong length")

  if(is.numeric(value)) {
    # entries in 'value' are integers
    if(!all(value %in% seq(TL)))
      stop("Tank numbers do not match tanklist")
    DF$tankid <- value
  } else {
    if(!any(nzchar(names(TL))))
      stop(paste("The tanklist does not assign names to the tanks;",
                 "tanks must be indexed by numbers"))
    # entries in 'value' are character or factor;
    # tanklist has names.
    if(is.character(value)) {
      if(!all(value %in% names(TL)))
        stop("Tank names do not match names in tanklist")
      DF$tankid <- factor(value, levels=names(TL))
    } else if(is.factor(value)) {
      if(!all(as.character(levels(value)) %in% names(TL)))
        stop("Tank names do not match names in tanklist")
      DF$tankid <- value
    } else stop("Unrecognised format for replacement value")
  }
  if(checkcollapse && length(unique(DF$tankid)) == 1)
    warning("Entire dive is conducted on a single tank")
  # update gas fractions
  DF <- reconcile.df(DF, TL)
  d$data <- DF
  return(d)
}
  
##################################################################
# methods for print, plot etc


plot.dive <- function(x, ...,
                      main=deparse(substitute(x)),
                      key.gases=c("text", "legend", "none"),
                      text.cex=1,
                      text.verticals=TRUE,
                      col.gases=1:length(tanklist(x)),
                      legendpos=c("top","bottom", "left", "right",
                        "topleft", "topright", "bottomleft", "bottomright",
                        "center"))
 {
  legendpos <- match.arg(legendpos)
  times <- times.dive(x)
  depths <- depths.dive(x)
  if(mean(diff(times)) <= 1 && missing(key.gases)) {
    key.gases <- "legend"
  } else key.gases <- match.arg(key.gases)

  # create plot 
  add <- resolve.defaults(list(...), list(add=FALSE))$add
  if(!add)
    do.call("plot.default",
            resolve.defaults(list(times, -depths, type="n"),
                             list(...),
                             list(main=main, axes=FALSE,
                                  xlab="Time (minutes)",
                                  ylab="Depth (metres)", lwd=2)))
  # plot axes
  axis(1, at=pretty(range(times)))
  yp <- pretty(range(depths))
  axis(2, at=-yp, labels=paste(yp))
  
  # determine plot colour for each gas
  ntank <- length(x$tanklist)
  conformit <- function(argu, n, default=1) {
    arguname <- deparse(substitute(argu))
    if(is.null(argu))
      argu <- default
    if(length(argu) == n)
      return(argu)
    if(length(argu) == 1)
      return(rep(argu, n))
    stop(paste(arguname, "has the wrong length"))
  }
  col.gases <- conformit(col.gases, ntank)
  
  # plot segments  
  noway <- nrow(x$data)
  tankid <- as.integer(x$data$tankid)
  do.call("segments",
          resolve.defaults(list(times[-noway], -depths[-noway],
                                times[-1], -depths[-1]),
                           list(...),
                           list(col=col.gases[tankid[-noway]],
                                lwd=2)))

  # annotate
  if(!x$airdive) 
    switch(key.gases,
           none = {},
           legend = {
             ug <- summary(x)$uniquegases
             legend(legendpos, lty=1, col=col.gases, legend=ug)
           },
           text = {
             labels <- with(x$data, gasnames(fO2, fN2, fHe))[-nrow(x$data)]
             rates <- abs(diff(depths)/diff(times))
             aspect <- diff(range(depths))/diff(range(times))
             vertical <- (rates > 2 * aspect)
             n <- length(times)
             midx <- (times[-n] + times[-1])/2
             midy <- -x$profile(midx)
             text(midx[!vertical],  midy[!vertical], labels[!vertical],
                  pos=3, cex=text.cex)
             if(text.verticals)
               text(midx[vertical], midy[vertical], labels[vertical],
                    srt=90, pos=4, offset=1, cex=text.cex)
           })
  invisible(NULL)
}

print.dive <- function(x, ..., seconds=TRUE) {
  cat("Dive profile\n")
  depths <- depths.dive(x)
  times <- times.dive(x)
  if(seconds) {
    mins <- floor(times)
    secs <- floor((times * 60) %% 60)
    watchtimes <- paste(mins, ":",
                        ifelse(secs < 10, "0", ""),
                        secs, sep="")
  } else 
    watchtimes <- paste(round(times))
  if(any(iii <- is.infinite(times)))
    watchtimes[iii] <- "Inf"
  
  if(x$airdive) {
    # air dive: print time and depth profiles
    cat("gas: air\n")
    print(data.frame(time=watchtimes, depth=depths))
  } else if(x$tanks.given) {
    # dive using specified tanks
    print(data.frame(time=watchtimes, depth=depths, tank=factor(x$data$tankid),
                     stringsAsFactors=FALSE), quote=FALSE)
    cat("\nTank contents:\n")
    tk <- x$tanklist
    gasblurbs <- unlist(lapply(tk, as.character, full=TRUE))
    tanknames <- names(tk)
    if(!any(nzchar(tanknames))) tanknames <- paste(1:length(tk))
    for(i in 1:length(gasblurbs))
      cat(paste(tanknames[i], ":\t", gasblurbs[i], "\n", sep=""))
  } else {
    # gases specified along the way.
    gases <- with(x$data, gasnames(fO2, fN2, fHe))
    if(length(unique(gases)) == 1) {
      cat(paste("gas:", gases[1], "\n"))
      print(data.frame(time=watchtimes, depth=depths))
    } else
    print(data.frame(time=watchtimes,
                     depth=depths,
                     gas=factor(gases),
                     stringsAsFactors=FALSE), quote=FALSE)
  }
  return(invisible(NULL))
}

summary.dive <- function(object, ...) {
  depths <- depths.dive(object)
  times <- times.dive(object)
  n <- length(depths)
  totaltime <- times[n]
  durations <- diff(times)
  middepths <- (depths[-1] + depths[-n])/2
  meandepth <- sum(durations * middepths)/totaltime
  gases <- with(object$data, gasnames(fO2, fN2, fHe))
  flat <- (diff(depths) == 0 & durations >= 1)
  stages <- data.frame(depth=(depths[-n])[flat],
                       time=durations[flat],
                       gas=(gases[-n])[flat])
  ugases <-
    if(object$airdive) "air" else
    if(!object$tanks.given) unique(gases) else 
    paste(names(object$tanklist),
          ":",
          unlist(lapply(object$tanklist, as.character)))

  # could be a double dive, etc
  ndives <- with(stages, 1 + sum(depth == 0))
  
  z <- list(depths=depths,
            times=times,
            totaltime=totaltime,
            stages=stages,
            ndives=ndives,
            maxdepth=max(depths),
            meandepth=meandepth,
            airdive=object$airdive,
            gasnames=gases,
            uniquegases=ugases)
  class(z) <- c("summary.dive", class(z))
  return(z)
}
  
print.summary.dive <- function(x, ...) {
  descrip <- switch(x$ndives,
                    '1' = "Dive to",
                    '2' = "Double dive to maximum depth",
                    paste("Sequence of", x$ndives, "dives to maximum depth"))
  cat(paste(descrip, x$maxdepth, "metres"))
  if(x$airdive)
    cat(" on air\n")
  else {
    gases <- unique(x$gasnames)
    ngas <- length(gases)
    if(ngas <= 5) 
      cat(paste("\n", ngettext(ngas, "Gas: ", "Gases: "),
                paste(gases, collapse=", "), "\n", sep=""))
    else
      cat("More than 5 gas settings")
  }
  cat(paste("Total dive time:", round(x$totaltime,1), "minutes\n"))
  if(nrow(x$stages) > 0) {
    cat("Stages:\n")
    if(x$airdive)
      print(x$stages[,1:2])
    else
      print(x$stages)
  }
  cat(paste("Mean depth", round(x$meandepth,1), "metres\n"))
  return(invisible(NULL))
}
  
###################################################################
#  Ascent/descent rates/times
#

is.rate <- function(x) { inherits(x, "rate") }

rate <- function(speed=NULL, time=NULL, up=TRUE) {
  ngiven <- (!is.null(speed)) + !is.null(time)
  if(ngiven != 1)
    stop("Exactly one of \"speed\" and \"time\" should be specified")
  if(!is.null(speed) && speed <= 0)
    stop("speed should be a positive number (metres/minute)")
  if(!is.null(time) && time <= 0)
    stop("time should be a positive number (minutes)")
  ra <- list(speed=speed, time=time, up=up)
  class(ra) <- c("rate", class(ra))
  return(ra)
}

ascent <- function(speed=NULL, time=NULL) { return(rate(speed, time, TRUE)) }

descent <- function(speed=NULL, time=NULL) { return(rate(speed, time, FALSE)) }

print.rate <- function(x, ...) {
  cat(paste(
            if(x$up) "ascent" else "descent",
            if(!is.null(x$speed)) 
              paste("rate", x$speed, "metres/minute")
            else
              paste("time", x$time, "minutes"),
            "\n"))
  return(invisible(NULL))
}

timetaken <- function(start, finish, uprate, downrate) {
  if(!is.rate(uprate) || !uprate$up)
    stop("\"uprate\" should be an ascent rate")
  if(!is.rate(downrate) || downrate$up)
    stop("\"downrate\" should be a descent rate")
  rat <- if(start > finish) uprate else downrate
  if(!is.null(rat$time))
    return(rat$time)
  else
    return(abs(finish-start)/rat$speed)
}


#################################################################
#
# Manipulating dives
#

chop.dive <- function(d, t0=0, t1=max(times.dive(d))) {
  stopifnot(is.dive(d))
  stopifnot(t1 > t0)
  times <- times.dive(d)
  df <- d$data
  if(t0 < 0)
    stop("The start time t0 is before the beginning of the dive")
  if(t1 > max(times))
    stop("The termination time t1 is later than the end of the dive")

  intervening <- (times >= t0 & times <= t1)
  middle <- any(intervening)
  if(middle) {
    # extract data from waypoints between t0 and t1 inclusive
    ra <- range(which(intervening))
    first <- ra[1]
    last <- ra[2]
    dfOUT <- df[first:last, ]
  } else {
    # there are no intervening waypoints
    dfOUT <- NULL
    # index for data that apply to dive
    first <- last <- max(which(times <= t0))
  }
  
  if(times[first] != t0) {
    # start time is between two waypoints
    # interpolate to find depth
    dep0 <- d$profile(t0)
    # tack on an initial row
    df0 <- df[first,]
    df0$time <- t0
    df0$depth <- dep0
    dfOUT <- rbind(df0, dfOUT)
  }
  if(times[last] != t1) {
    # finish time is between two waypoints
    # interpolate to find depth
    dep1 <- d$profile(t1)
    # tack on a final row
    df1 <- df[last,]
    df1$time <- t1
    df1$depth <- dep1
    dfOUT <- rbind(dfOUT, df1)
  }
  # reset the clock
  dfOUT$time <- with(dfOUT, time - time[1])
  # create dive profile 
  pro <- approxfun(dfOUT$time, dfOUT$depth,
                   method="linear", rule=2, ties="ordered")
  airdive <- all(with(dfOUT, gasnames(fO2, fN2, fHe)) == "air")
  result <- d
  result$profile <- pro
  result$data    <- dfOUT
  result$airdive <- airdive
  return(result)
}

dive.segment <- function(d, i) {
  # make a `dive' out of one segment (between two waypoints) of a dive d
  df <- d$data[c(i, i+1), ]
  airdive <- all(with(df, gasnames(fO2, fN2, fHe)) == "air")
  pro <- approxfun(df$time, df$depth, method="linear", rule=2, ties="ordered")
  result <- d
  result$profile <- pro
  result$data    <- df
  result$airdive <- airdive
  return(result)
}


###########################################################
# tank lists

tanklist <- function(d) {
  stopifnot(is.dive(d))
  return(d$tanklist)
}

"tanklist<-" <- function(d, value) {
  if(!is.dive(d))
    stop("In tanklist(d) <- value, d must be a dive object")

  if(is.null(value)) {
    # wipe all tank data
    d$tanks.given <- FALSE
    d$tanklist <- list(air)
    d$airdive <- TRUE
    d$data$tankid <- 1
    d$data$fO2 <- air$fO2
    d$data$fN2 <- air$fN2
    d$data$fHe <- air$fHe
    return(d)
  }
  if(!is.list(value) ||
     length(value) == 0 ||
     !all(unlist(lapply(value, is.gas))))
    stop("In tanklist(d) <- value, value must be a list of gases")

  tk <- d$tanklist
  if(is.null(tk)) {
    # Dive does not currently contain any tank information
    # Assign tanklist
    d$tanklist <- value 
    # Assume tank 1 used throughout dive
    tanknames <- names(value)
    if(all(nzchar(tanknames)))
      d$data$tankid <- factor(tanknames[1], levels=tanknames)
    else 
      d$data$tankid <- 1
    # Fill out gas fractions
    g <- tanklist[[1]]
    d$data$fO2 <- g$fO2
    d$data$fN2 <- g$fN2
    d$data$fHe <- g$fHe
    # Was it an air dive?
    d$airdive <- is.air(g)
    return(d)
  }
  
  if(length(tk) != length(value)) {
    # check whether existing tankid values would be valid
    mold <- max(as.integer(d$data$tankid))
    mnew <- length(value)
    if(mold > mnew)
      stop(paste("New tanklist contains only", mnew, "tanks;",
                 "dive schedule requires", mold, "tanks"))
  }
  
  # reconcile names of tanks
  hasoldnames <- any(nzchar(names(tk)))
  hasnewnames <- any(nzchar(names(value)))
  if(!hasnewnames && hasoldnames)
    names(value) <- names(tk)

  d$tanks.given <- TRUE
  d$tanklist <- value
  df <- d$data
  # update gas fractions
  df <- reconcile.df(df, value)

  if(hasnewnames) 
    df$tankid <- factor(names(value)[as.integer(df$tankid)],
                        levels=names(value))

  d$data <- df
  d$airdive <- all(unlist(lapply(value, is.air)))
  return(d)
}

allspecies <- function(d, inert=TRUE) {
  stopifnot(is.dive(d))
  gf <- d$data[, c("fO2", "fN2", "fHe")]
  present <- (apply(as.matrix(gf), 2, max) > 0)
  spec <- c("O2", "N2", "He")
  out <- spec[as.logical(present)]
  if(inert)
    out <- out[out != "O2"]
  return(out)
}

reconcile.df <- function(df, tanks) {
  ## Ensure gas fractions are correctly determined by tank id
  tkid <- df$tankid
  ntanks <- length(tanks)
  if(is.numeric(tkid)) {
    if(!all(unique(tkid) %in% seq_len(ntanks)))
      stop("Tank numbers out-of-bounds")
    for(k in seq_len(ntanks)) {
      relevant <- (tkid == k)
      if(any(relevant)) {
        g <- tanks[[k]]
        df$fO2[relevant] <- g$fO2
        df$fN2[relevant] <- g$fN2
        df$fHe[relevant] <- g$fHe
      }
    }
  } else {
    tkname <- as.factor(tkid)
    nama <- names(tanks)
    if(!all(levels(tkname) %in% nama))
      stop("Tank names not recognised")
    for(k in seq_len(ntanks)) {
      relevant <- (tkname == nama[k])
      if(any(relevant)) {
        g <- tanks[[k]]
        df$fO2[relevant] <- g$fO2
        df$fN2[relevant] <- g$fN2
        df$fHe[relevant] <- g$fHe
      }
    }
  }
  return(df)
}

# Resolve conflicts between several sets of defaults
# Usage:
#     resolve.defaults(list1, list2, list3, .......)
# where the earlier lists have priority 
#
resolve.defaults <- function(...) {
  arglist <- list(...)
  argue <- list()
  if((n <- length(arglist)) > 0)  {
    for(i in seq(n))
      argue <- append(argue, arglist[[i]])
  }
  if(!is.null(nam <- names(argue))) {
    named <- (nam != "")
    arg.unnamed <- argue[!named]
    arg.named <-   argue[named]
    if(any(discard <- duplicated(names(arg.named)))) 
      arg.named <- arg.named[!discard]
    argue <- append(arg.unnamed, arg.named)
  }
  return(argue)
}
