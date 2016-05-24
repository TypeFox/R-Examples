summary.mult.ev.data.sim <-
function(object, ...)
  {
    if(!inherits(object, "mult.ev.data.sim")) stop("Wrong data type")
    sub.risk    <- vector()
    num.events  <- vector()
    foltime     <- vector()
    med.foltime <- vector()
    mean.ep.sub <- vector()
    dens.incid  <- vector()
    ev.num      <- vector()
    for (i in 1:attr(object,"nsit"))
    {
       ev.num[i]      <- i
       sub.risk[i]    <- attr(object,"n")
       num.events[i]  <- as.integer(sum(object$status[object$ev.num==i]))
       foltime[i]     <- sum(object$time[object$ev.num==i])
       med.foltime[i] <- median(object$time[object$ev.num==i]) 
       dens.incid[i]  <- num.events[i]/foltime[i]
    }
    ans <- data.frame(ev.num, sub.risk, num.events, 
                 foltime, med.foltime, dens.incid)
    class(ans)  <- "summary.mult.ev.data.sim"
    return(ans)
  }
