### data

mapfit.time.data.frame <- function(time, difftime) {
  if (missing(time)) {
    if (missing(difftime)) {
      stop("error: mapfit.time.data.frame")
    }
  } else {
    time <- sort(time)
    if (time[1] != 0) {
      difftime <- diff(c(0, time))
    } else {
      difftime <- diff(time)
    }
  }
  new("mapdata.time", size=length(difftime),
    data=data.frame(time = difftime, counts = 0, instant = 1))
}

mapfit.group.data.frame <- function(counts, breaks, difftime, instant) {
  if (missing(instant)) {
    instant <- array(0, length(counts))
  }
  if (missing(breaks)) {
    if (missing(difftime)) {
      breaks <- 0:length(counts)
    } else {
      breaks <- c(0,cumsum(difftime))
    }
  }
  if (breaks[1] != 0) {
    breaks <- c(0, breaks)
    counts <- c(NA, counts)
    instant <- c(0, instant)
  }
  new("mapdata.group", size=length(counts),
    data=data.frame(time = diff(breaks), counts = counts, instant = instant))
}

setAs("mapdata.time", "mapdata.group", function(from, to) {
  new("mapdata.group", size=from@size, data=from@data)
  })

  ###########

##setMethod("[", "mapdata", function(x, i, j, ..., drop) callGeneric(x@data, i, j, ..., drop))
setMethod("print", signature(x = "mapdata"), function(x, ...) print(x@data))
setMethod("summary", signature(object = "mapdata"), function(object, ...) summary(object@data))

###########

setMethod("mean", signature(x = "mapdata"), function(x, ...) mapfit.mean(x))

setMethod("mapfit.mean", signature(x = "mapdata"),
  function(x, ...) {
    dat <- x@data
    dat <- dat[is.finite(dat$time) & is.finite(dat$counts) & is.finite(dat$instant),]
    time <- cumsum(dat$time)
    ttime <- time[dat$instant == 1]
    m <- sum(time * dat$counts) + sum(ttime)
    m / (sum(dat$counts) + length(ttime))
  }
)

