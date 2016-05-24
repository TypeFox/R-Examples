calc.v <- function(event, time) {
  ## make event a logical vector
  event <- as.logical(event)
  
  ## Set NA's to 0 to stop buggyness
  time[is.na(time)] <- 0L
  event[is.na(event)] <- FALSE

  orderTime <- order(time)
  sortTime <- time[orderTime]
  sortEvent <- event[orderTime]

  ## determin number of events
  n.select <- length(event)
  
  ## data.frame elements were
  ##  events$time : junk
  ##  events$rank : H
  ##  events$diff.no : dH0
  ##  events$diff.yes : dH1
  ##  events$rank.adj : Hadj
  events <- data.frame(time = unique(sortTime),
                       rank = unique(rank(sortTime,
                         ties.method="max")/n.select))
  events$diff.neg <- I(0.0)
  events$diff.pos <- I(0.0)

  ## get rid of all div by zero errors.
  ## this is okay because where ever rank is 1 the numerator will also be 0
  ## and then the value of that should be zero
  events$rank.adj <- ifelse(events$rank == 1L, 0.99, events$rank)

  ## find time values for patients that are selected but did not have
  ## an event
  ##
  ## data.frame elements were
  ##  events.neg$time : junk0
  ##  events.neg$diff : dH0a
  temp.var <- sortTime[!sortEvent]
  events.neg <- data.frame(time = unique(temp.var),
                          diff = diff(c(0L, unique(rank(temp.var,
                            ties.method="max")/n.select))))
  
  ## find time values for patients that are selected but did have an event
  ##
  ## data.frame elements were
  ##  events.pos$time : junk1
  ##  events.pos$diff : dH1a
  temp.var <- sortTime[sortEvent]
  events.pos <- data.frame(time = unique(temp.var),
                           diff = diff(c(0L, unique(rank(temp.var,
                             ties.method="max")/n.select))))

  # expand the diffs for event and no event to the full events$time vector length
  events$diff.neg[match(events.neg$time, events$time)] <- events.neg$diff
  events$diff.pos[match(events.pos$time, events$time)] <- events.pos$diff

  
  gamma0 <- exp(sapply(events$time, data.var=events,
                       FUN=function(time, data.var) {
                         sum((data.var$time < time) * data.var$diff.neg /
                             (1L - data.var$rank.adj))
                       }))
  
  calc.vrow <- function(timeval, data.var, gamma0, time, event) {
    phi <- data.var$time <= timeval

    gamma1 <- sapply(seq_len(nrow(data.var)), data.var=data.var, gamma0=gamma0,
                     phi=phi,
                     FUN=function(i, data.var, gamma0, phi) {
                       sum((data.var$time[i] < data.var$time) * phi *
                           gamma0 * data.var$diff.pos) /
                             (1L - data.var$rank.adj[i])
                     })

    integralc <- .C("integral", as.integer(nrow(data.var)),
                    as.double(data.var$time), as.double(data.var$time),
                    as.double(data.var$diff.neg), as.double(data.var$diff.pos),
                    as.double(data.var$rank.adj), as.double(gamma0),
                    as.double(phi), gamma2=double(nrow(data.var)), PACKAGE='sensitivityPStrat')

    gamma2 <- integralc[["gamma2"]]

    index <- match(time, data.var$time, nomatch=0L)
    
    v <- phi[index] * gamma0[index] * event + gamma1[index] * (!event) - gamma2[index]

    return(v)
  }
  
  v <- t(sapply(events.pos$time, data.var=events, gamma0=gamma0, time=time,
                event=event, FUN=calc.vrow))

  return(v)
}
