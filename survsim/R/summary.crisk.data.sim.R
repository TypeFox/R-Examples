summary.crisk.data.sim <-
function(object, ...)
  {
    if(!inherits(object, "crisk.data.sim")) stop("Wrong data type")
    sub.risk    <- vector()
    num.events  <- vector()
    foltime     <- vector()
    dens.incid  <- vector()
    cause       <- vector()
    for (i in 1:attr(object,"nsit"))
    {
       cause[i]       <- i
       sub.risk[i]    <- attr(object,"n")
       num.events[i]  <- as.integer(sum(!is.na(object$status[object$cause==i])))
       foltime[i]     <- sum(object$time)
       dens.incid[i]  <- num.events[i]/foltime[i]
    }
    ans <- data.frame(cause, sub.risk, num.events, dens.incid)
    class(ans)  <- "summary.crisk.data.sim"
    return(ans)
  }
