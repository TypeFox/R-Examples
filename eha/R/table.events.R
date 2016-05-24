table.events <- function(enter = rep(0, length(exit)),
                         exit,
                         event,
                         strict = TRUE)
{
  n <- length(exit)

  ## Check input data:
  if ( length(enter) != n ) stop("enter and exit must have equal length.")
  if ( length(event) != n ) stop("event and exit must have equal length.")
  ##
  
  event <- (event != 0) ## 0 (FALSE) = censoring, else (TRUE) = event

  times <- c(unique(sort(exit[event])))
  nn <- length(times)

  rs.size <- double(nn)
  n.events <- double(nn)

  for (i in 1:nn) ## Try to avoid this loop!
    {
      rs.size[i] <- sum( (enter < times[i]) &
                        (exit >= times[i]) )
      n.events[i] <- sum( (exit == times[i]) & event )
    }

  stop.at <- which(rs.size == n.events)
  if (strict & length(stop.at))
    {
      stop.at <- min(stop.at) - 1
      if (stop.at <= 0)
          stop("First risk set is all events! Try 'strict = FALSE'")
      times <- times[1:stop.at]
      n.events <- n.events[1:stop.at]
      rs.size <- rs.size[1:stop.at]
    }
      
  return ( list(times         = times,
                events        = n.events,
                riskset.sizes = rs.size)
          )
}

