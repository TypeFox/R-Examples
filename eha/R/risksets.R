risksets <- function (x, strata = NULL, max.survs = NULL, members = TRUE){
    ## x is a Surv (survival) object.

  nn <- NROW(x)
  if (is.null(strata)){
      strata <- rep(1, nn)
  }else{
      if (length(strata) != nn) stop("'strata' has wrong length")
      else
          strata <- as.integer(factor(strata))
  }

  if (is.null(max.survs)) max.survs <- nn - 1
  if (NCOL(x) == 2){
      enter <- numeric(nn)
      exit <- x[, 1]
      event <- (x[, 2] != 0)
  }else{
      if (NCOL(x) != 3) stop("'x' is not a Surv object")
      enter <- x[, 1]
      exit <- x[, 2]
      event <- (x[, 3] != 0) ## TRUE == event
  }

  ord <- order(strata, exit, -event)
  strata <- strata[ord]
  enter <- enter[ord]
  exit <- exit[ord]
  event <- event[ord]
  w.totrs <- sum(nn) ## Working 'totrs'
  ns <- max(strata)
  nstra <- c(0, cumsum(table(strata)))
  
  counts <- .C("sizes",
               as.integer(ns),
               as.integer(nn), 
               as.double(enter),
               as.double(exit),
               as.integer(event),
               ##
               antrs = integer(ns),
               as.integer(nstra),
               risktimes = double(w.totrs),
               ##
               n.events = integer(w.totrs),
               size = integer(w.totrs),
               totrs = integer(1),
               ## DUP = FALSE,
               PACKAGE = "eha"
               )

  counts$risktimes <- counts$risktimes[1:counts$totrs]
  counts$n.events <- counts$n.events[1:counts$totrs]
  counts$size <- counts$size[1:counts$totrs]

  totsize <- sum(counts$size)
  if (totsize >= 2^31) stop("Too large  risk sets.") 
  totevents <- sum(counts$n.events)

  Eventset <- NULL
  Riskset <- NULL
  sample_fraction <- NULL
  
  if (members){
      res <- .C("risk_get",
                as.integer(max.survs),
                as.integer(nn),
                as.integer(ns),
                ##
                as.double(enter),
                as.double(exit),
                as.integer(event),
                ##
                as.integer(nstra),
                as.integer(length(nstra)),
                ##
                new.totrs = integer(1),  ## If sampling...
                ##
                as.integer(counts$antrs),
                as.integer(counts$n.events),
                size = as.integer(counts$size), ## If sampling...
                as.double(counts$risktimes),
                eventset = integer(totevents),
                riskset = integer(totsize),
                ## DUP = FALSE,
                PACKAGE = "eha"
                )
      Size <- res$size ## Previously out-commented; why??
      N <- counts$size - counts$n.events
      sample_fraction <- numeric(length(Size))
      sample_fraction[N > 0] <- (Size - counts$n.events) /
          (counts$size - counts$n.events)
      sample_fraction[N <= 0] <- 1 # No survivors!
      Eventset <- ord[res$eventset]
      Riskset <- ord[res$riskset[1:res$new.totrs]]
  }

  rs <- list(ns = ns,
             antrs = counts$antrs,
             risktimes = counts$risktimes,
             n.events = counts$n.events,
             ##size = counts$size,
             size = Size,
             eventset = Eventset,
             riskset = Riskset,
             sample_fraction = sample_fraction)
  class(rs) <- "risksets"
  
  invisible(rs)
}
