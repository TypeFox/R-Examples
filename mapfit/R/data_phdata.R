
phfit.time.data.frame <- function(time, weights) {
  if (missing(weights)) {
    weights <- array(1, length(time))
  }
  # d <- data.frame(time = diff(c(0, sort(time))),
  #                    weights = weights[order(time)])
  tt <- sort(time)
  new("phdata.wtime", size=length(time),
    data=data.frame(time = tt, weight = weights[order(time)]),
    diff=diff(c(0, tt)))
}

phfit.group.data.frame <- function(counts, breaks, difftime, instant) {
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
  new("phdata.group", size=length(counts),
    data=data.frame(time = diff(breaks), counts = counts, instant = instant))
}

##setMethod("[", "phdata", function(x, i, j, ..., drop) callGeneric(x@data, i, j, ..., drop))
setMethod("print", signature(x = "phdata"), function(x, ...) print(x@data))
setMethod("summary", signature(object = "phdata"), function(object, ...) summary(object@data))

###########

setMethod("mean", signature(x = "phdata"), function(x, ...) mapfit.mean(x))

setMethod("mapfit.mean", signature(x = "phdata.wtime"),
  function(x, ...) {
    dat <- x@data
    sum(dat$time * dat$weight) / sum(dat$weight)
  }
)

setMethod("mapfit.mean", signature(x = "phdata.group"),
  function(x, ...) {
    dat <- x@data
    dat <- dat[is.finite(dat$time) & is.finite(dat$counts) & is.finite(dat$instant),]
    time <- cumsum(dat$time)
    ttime <- time[dat$instant == 1]
    m <- sum(time * dat$counts) + sum(ttime)
    m / (sum(dat$counts) + length(ttime))
  }
)

