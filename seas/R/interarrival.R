"interarrival" <-
function(x, var="precip", p.cut=0.3, inv=FALSE) {
  orig <- as.character(substitute(x))[[1]]
  sc <- seas.df.check(x, orig, var)
  x <- x[!is.na(x[,var]),]
  x$diff <- c(0, as.integer(diff(x$date)))
  x$wet <- (x[,var] > p.cut)
  if (inv)
    x$wet <- !x$wet
  j <- 0
  date <- x$date[length=0]
  wet <- dry <- integer(0)
  was.wet <- NA
  for (i in 1:nrow(x)) {
    if (x$diff[i] != 1) {
      mode <- TRUE
      if (!is.na(was.wet)) {
        if (was.wet) {
          wet[j] <- NA
        } else {
          dry[j] <- NA
        }
      }
    } else {
      if (mode) {
        j <- j + 1
        mode <- FALSE
      }
      if (x$wet[i]) {
        if (!was.wet) {
          date[j] <- x$date[i]
          wet[j] <- 0
        }
        wet[j] <- wet[j] + 1
      } else {
        if (was.wet) {
          j <- j + 1
          date[j] <- NA
          dry[j] <- 0
          wet[j] <- NA
        }
        dry[j] <- dry[j] + 1
      }
    }
    was.wet <- x$wet[i]
  }
  if (was.wet)
    wet[j] <- NA
  else
    dry[j] <- NA
  if (inv) {
    inv <- wet
    wet <- dry
    dry <- inv
  }
  s <- (!is.na(date) & !(is.na(wet) & is.na(dry)))
  d <- data.frame(date=date[s], dry=dry[s], wet=wet[s])
  class(d) <- c("interarrival", "data.frame")
  attr(d, "id") <- sc$id
  attr(d, "name") <- sc$name
  attr(d, "year.range") <- sc$year.range
  attr(d, "p.cut") <- p.cut
  attr(d, "inv") <- inv
  d
}
