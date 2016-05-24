library(zoo)

# Upon input
#   z is a zoo object containing input data. E.g. this could be all the 
#     prices of a bunch of stocks. The column name is the unit name.
#   events is a data.frame containing 2 columns. The first column
#     ("unit") is the name of the unit. The second column is the date/time
#     ("when") when the event happened.
# For each event, the outcome can be:
#   unitmissing : a unit named in events isn't in z
#   wrongspan : the event date isn't placed within the span of data for the unit
#   wdatamissing: too many NAs within the crucial event window.
#   success : all is well.
# A vector of these outcomes is returned.
phys2eventtime <- function(z, events, width=10) {
  # Just in case events$unit has been sent in as a factor --
  events$unit <- as.character(events$unit)
  if(is.factor(events$when)) stop("Sorry you provided a factor as an index")
  # Given a zoo time-series vector x, and an event date "when",
  # try to shift this vector into event time, where the event date
  # becomes 0 and all other dates shift correspondingly.
  # If this can't be done, then send back NULL with an error code.
  timeshift <- function(x, when) {
    location <- findInterval(when, index(x))
    if ((location <= 1) | (location >= length(x))) {
      return(list(result=NULL, outcome="wrongspan"))
    }
    remapped <- zoo(as.numeric(x), order.by=(-location+1):(length(x)-location))
    list(result=remapped, outcome="success")
  }

  # Main loop to build up a data object in event time --
  outcomes <- character(nrow(events))
  z.e <- zoo(1, order.by=as.integer(1)) # zoo::cbind() requires initialising z.e
  for (eventnum in 1:nrow(events)) {
    if (!(events$unit[eventnum] %in% colnames(z))) {
      outcomes[eventnum] <- "unitmissing"
      next
    }
    attempt <- timeshift(z[,events$unit[eventnum]], events$when[eventnum])
    if (attempt$outcome=="success") {
      z.e <- cbind(z.e, attempt$result)
    }
    outcomes[eventnum] <- attempt$outcome
  }
  outcomes <- outcomes
  z.e <- z.e[,-1, drop = FALSE]                       #get rid of that junk initialisation
  colnames(z.e) <- which(outcomes=="success")
  ## Now worry about whether there's information within the event window
  ## (This entire cleaning+checking can be switched off by specifying width=0)
  badcolumns <- NULL
  if (width > 0) {
    for (i in 1:ncol(z.e)) {
      tmp <- z.e[,i]
      tmp <- na.locf(tmp, na.rm=FALSE, maxgap=4)
      tmp <- na.locf(tmp, na.rm=FALSE, maxgap=4, fromLast=TRUE)
      tmp2 <- window(tmp, start=-width, end=+width)
      if (any(is.na(tmp2))) {
        outcomes[as.numeric(colnames(z.e)[i])] <- "wdatamissing"
        badcolumns <- c(badcolumns, i)
      } else {
        z.e[,i] <- tmp                # Put the fixed up column back in.
      }
    }
    if (any(outcomes == "wdatamissing")) {
      z.e <- z.e[, -badcolumns]
    }
  }
  # Check that we're okay
  stopifnot(sum(outcomes=="success") == NCOL(z.e))
  list(z.e=z.e, outcomes=factor(outcomes))
}
