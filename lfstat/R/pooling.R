
find_droughts <- function(x, threshold = vary_threshold, ...) {
  if(!inherits(x, "xts")) x <- as.xts(x)
  x <- .check_xts(x)

  discharge <- if(ncol(x) == 1) x[, 1] else {
    if(!"discharge" %in% names(x)) {
      stop("'x' must eihter be an object of class lfobj, an xts object with one column or with a column named 'discharge'.")
    }
    x[, "discharge"]
  }

  att <- xtsAttributes(x)
  # we want volumes always in qubic metre and times always in seconds
  f.vol <- .conv_factor(from = att$unit.parsed["volume"], to = "m",
                        dimension = "vol")
  f.time <- att$deltat * .conv_factor(from = att$unit.parsed["time"],
                                      to = "secs", dimension = "time")
  f <- f.time / f.vol

  names(discharge) <- "discharge"

  if(is.function(threshold)) threshold <- threshold(discharge, ...)

  start <- 1
  event.no <- 0
  len <- nrow(x)

  x <- cbind(discharge = discharge,
             threshold  = threshold,
             def.increase = as.vector(threshold - coredata(discharge)) * f,
             event.no = numeric(len))

  deficit <- x$def.increase >= 0
  x$event.no <- group(deficit, new.group.na = TRUE)

  x$event.no[deficit] <- as.numeric(factor(x$event.no[deficit]))
  x$event.no[!deficit | is.na(x$def.increase)] <- 0

  class(x) <- unique(c("deficit", class(x)))
  return(x)
}

pool_ic <- function(x, tmin = 5, ratio = 0.1) {
  tab <- summary(x, drop_minor = 0)

  x$event.orig <- x$event.no
  if(nrow(tab)) {

    # inter event time and volume
    ti <- as.numeric(difftime(tail(tab$start, -1), head(tab$end, -1),
                              units = "days"))

    vi <- numeric(length(nrow(tab) - 1))
    for(i in seq_len(nrow(tab) - 1)) {
      start <- tab$end  [i]   + as.difftime(1, units = "days")
      end   <- tab$start[i+1] - as.difftime(1, units = "days")
      vi[i] <- sum(x[paste(start, end, sep = "/"), "def.increase"], na.rm = T)
    }

    tab$vol.pooled <- tab$volume
    event <- c(1, 2)

    repeat {
      if(ti[event[2] - 1] < tmin && abs(vi[event[2] - 1] / tab$vol.pooled[event[2] - 1]) < ratio){
        tab$event.no[event[2]] <- event[1]
        tab$vol.pooled[event[2]] <- sum(tab$vol.pooled[event[2] - 0:1]) + vi[event[2] - 1]

        rng <- with(tab, paste(end[event[2] - 1], end[event[2]], sep = "/"))
        x$event.no[rng] <- tab$event.no[event[1]]

        event[2] <- event[2] + 1
      } else {
        event <- event[2] + c(0, 1)
      }

      if(max(event) == nrow(tab)) break
    }
  }

  n.pooled <- length(unique(x$event.orig)) - length(unique(x$event.no))

  xtsAttributes(x) <- c(xtsAttributes(x), list(pooling = "IC",
                                               n.pooled = n.pooled))
  return(x)
}


pool_ma <- function(x, n = 10) {

  x$discharge <- ma(x$discharge, n = n, sides = 2)

  y <- find_droughts(x$discharge, threshold = x$threshold)
  y$event.orig <- x$event.no

  n.pooled <- length(unique(y$event.orig)) - length(unique(y$event.no))
  xtsAttributes(y) <- c(xtsAttributes(y), list(pooling = "Moving Average",
                                               n.pooled = n.pooled))

  return(y)
}

pool_sp <- function(x) {
  # ignores precomupted event numbers
  x$event.orig <- x$event.no
  x$event.no <- 0

  len <- nrow(x)
  start <- 1
  event.no <- 0
  repeat {
    # find start of drought, first day with a positive increase of the deficit
    start <- start + which(x$def.increase[start:len] > 0)[1] - 1
    if(is.na(start)) break

    # can't use head() because n = -0 returns an empty object
    deficit <- cumsum(x$def.increase[start:len])

    # range of the event
    # Tallaksen excluded the day of zero-crossing, maybe we keep it?
    end <- which(deficit < 0)[1] - 1
    if(is.na(end)) break
    rng <- seq(start, start + end - 1)

    event.no <- event.no + 1
    x$event.no[rng] <- event.no

    start <- start + end
  }

  n.pooled <- length(unique(x$event.orig)) - length(unique(x$event.no))
  xtsAttributes(x) <- c(xtsAttributes(x), list(pooling = "Sequent Peak",
                                               n.pooled = n.pooled))
  return(x)
}


summarize.drought <- function(x, drop_minor = c("volume" = 0, "duration" = 0),
                              duration = c("event", "peak")) {
  duration <- match.arg(duration)

  def.vol <- cumsum(as.numeric(x$def.increase))
  time.ind <- time(x)


  # for Sequent Peak algorithm
  # drought duration is defined as time until maximum depletion

  if (duration == "peak") {
    duration <- which.max(def.vol)
    time <- time.ind[duration]
  } else {
    duration <- length(time.ind)
    time <- time.ind[1]
  }

  # neglect minor events
  if (def.vol[duration] < drop_minor["volume"] ||
      duration < drop_minor["duration"] ) {
    return(data.frame())
  }

  data.frame(event.no = coredata(x$event.no)[1],
             start = time.ind[1],
             time = time,
             end = tail(time.ind, 1),
             volume = def.vol[duration],
             duration = duration)
}



.parse_minor_arg <- function(arg, x) {
  if(any(grepl("%", arg, fixed = T))) {
    x <- summary(x, drop_minor = c("volume" = 0, "duration" = 0))
  }

  f <- function(val, sample) {
    if(grepl("%", val, fixed = T)) {
      val <- as.numeric(sub("%", "", val, fixed = T))
      if (val < 0 || val > 100) stop("fraction must be between 0% and 100%.")
      y <- max(sample, na.rm = T) * val / 100
    } else {
      y <- as.numeric(val)
    }

    return(y)
  }

  if (length(arg) == 1){
    if (arg == 0) {
      arg <- c("volume" = 0, "duration" = 0)
    } else {
      stop("if of length one, argument 'drop_minor' has to be zero.")
    }
  } else if (length(arg) == 2){
    if(is.null(names(arg))) {
      names(arg) <- c("volume", "duration")
    } else {
      pos <- charmatch(names(arg), c("volume", "duration"))
      if (any(is.na(pos) | pos == 0)) {
        stop("names 'volume' and 'duration' not found in argument 'drop_minor'")
      } else {
        arg <- arg[pos]
        names(arg) <- c("volume", "duration")
      }
    }
  } else {
    stop("argument 'drop_minor' has to be of length one or two.")
  }

  res <- numeric(length(arg))
  names(res) <- names(arg)
  for(i in seq_along(arg)) {
    res[i] <- f(arg[i], sample = x[, names(arg[i])])
  }

  return(res)
}


summary.deficit <- function(object,
                            drop_minor = c("volume" = "0.5%", "duration" = 5),
                            ...) {

  drop_minor <- .parse_minor_arg(drop_minor, object)

  x <- object[object$event.no > 0, ]

  pooling <- xtsAttributes(object)$pooling
  duration <- if(!is.null(pooling) && pooling == "Sequent Peak") "peak" else "event"

  y <- lapply(split(x, x$event.no), summarize.drought, drop_minor = drop_minor,
              duration = duration)
  omitted <- sum(sapply(y, function(x) length(x) == 0))
  total <- length(y)

  y <- do.call(rbind, y)

  class(y) <- c("summaryDrought", class(y))
  attr(y, "deficit") <- c(xtsAttributes(x),
                          list(drop_minor = drop_minor, omitted = omitted,
                               total = total))
  return(y)
}

print.summaryDrought <- function(x, ...) {
  attlist <- attr(x, "deficit")
  cat("Summary of", if (!is.null(attlist$pooling)) "pooled", "droughts")

  if(!is.null(attlist)) {
    if (!is.null(attlist$river)) cat("\nRiver:", attlist$river)
    if (!is.null(attlist$station)) cat(" at", attlist$station)
    if (!is.null(attlist$pooling)) cat("\nPooling:", attlist$pooling)
    if (!is.null(attlist$n.pooled)) cat(", ", attlist$n.pooled, " were pooled", sep = "")

    cat("\nUnits: volumes in m\u00B3, duration in days")

    if(attlist$omitted) {
      with(attlist, cat("\n\nFiltered", omitted, "minor events of", total , "total.",
                        "\nOnly droughts with volume >=", drop_minor["volume"],
                        "m\u00B3 and duration >=", drop_minor["duration"], "days are",
                        "reported."))
    }
  }

  cat("\n\n")

  NextMethod(x)
}


print.deficit <- function(x, ...) {
  # pool oder event.no
  cat("streamflow defict", fill = TRUE)
  cat("data.frame with", nrow(x), "observations and",
      length(unique(x$event.no)), "events\n\n")

  NextMethod(x)
}

plot.deficit <- function(x, type = "dygraph", ...) {
  if (type == "dygraph") {
    plot.deficit_dygraph (x, ...)
  } else {
    plot.xts(x$discharge, type = type, ...)
    lines(x$threshold, col = 2, ...)
  }
}


plot.deficit_dygraph <- function(x, ...) {
  arg <- list(...)
  if("step" %in% names(arg)) step <- arg$step else step = TRUE

  is.drought  <- x$event.no != 0
  x$lwr <- x$upr <- x$discharge

  # expand xts object with cols upr and lwr, indicating deficits
  x$lwr[is.drought] <- with(x[is.drought], ifelse(threshold >= discharge, discharge, threshold))
  x$upr[is.drought] <- with(x[is.drought], ifelse(threshold <= discharge, discharge, threshold))


  # for some reaseon range of deficit is one element too large at each margin
  pos <- diff(coredata(x[, "event.no", T]) != 0)
  border <- logical(nrow(x))
  #   border[which(pos == 1) + 1] <- T
  border[which(pos == -1)] <- T
  x[border, c("lwr", "upr")] <- x$threshold[border]

  attlist <- xtsAttributes(x)
  title <- .char2html(with(attlist, paste("River", river, "at", station)))
  ylab <- .char2html(paste("Flow in", attlist$unit))
  p <- dygraph(x[, c("discharge", "lwr", "threshold", "upr")],
               main = title, ylab = ylab) %>%
    dyRangeSelector() %>%
    dySeries("discharge", stepPlot = step, drawPoints = TRUE, color = "darkblue") %>%
    dySeries(c("lwr", "threshold", "upr"), stepPlot = step, color = "red",
             strokePattern = "dashed")

  tbl <- summary(x, ...)

  if(length(tbl)) {
    ttip <- with(tbl,
                 paste0("volume: ", signif(volume, 3), "\n",
                        "duration: ", duration, " days\n",
                        "intensity: ", signif(volume/duration, 3)))


    for (i in seq_len(nrow(tbl))) {
      p <- dyShading(p, from = tbl$start[i], to = tbl$end[i], color = "lightgrey")
      p <- dyAnnotation(p, x = round(mean(c(tbl$start[i], tbl$end[i]))),
                        text = tbl$event.no[i],
                        tooltip = ttip[i], width = 30)
    }
  }
  return(p)
}

