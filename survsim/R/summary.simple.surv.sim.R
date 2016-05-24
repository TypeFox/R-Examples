summary.simple.surv.sim <-
function(object, ...)
  {
    if(!inherits(object, "simple.surv.sim")) stop("Wrong data type")
    sub.risk    <- vector()
    num.events  <- vector()
    foltime     <- vector()
    med.foltime <- vector()
    mean.ep.sub <- vector()
    dens.incid  <- vector()
    sub.risk[1]    <- attr(object,"n")
    num.events[1]  <- as.integer(sum(object$status))
    foltime[1]     <- sum(object$stop)
    med.foltime[1] <- median(object$stop) 
    mean.ep.sub[1] <- sum(object$status)/dim(object)[1]
    dens.incid[1]  <- num.events[1]/foltime[1]
    
    ans <- data.frame(sub.risk, num.events, mean.ep.sub, 
                      foltime, med.foltime,
                      dens.incid)
    class(ans)  <- "summary.simple.surv.sim"
    return(ans)
  }
