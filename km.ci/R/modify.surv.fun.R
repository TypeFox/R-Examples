"modify.surv.fun" <-
function(survi,start,end,method)
{
    # This function simply modifys an survival object
    # and cuts off the ends where no upper and lower
    # boundary is calculated
    survi <- survi
    survi$time <- survi$time[start:end]
    survi$n.risk <- survi$n.risk[start:end]
    survi$n.event <- survi$n.event[start:end]
    survi$surv <- survi$surv[start:end]
    survi$std.err <- survi$std.err[start:end]
    survi$upper <- survi$upper[start:end]
    #survi$upper[survi$upper>1] <- 1
    survi$lower <- survi$lower[start:end]
    #survi$lower[survi$lower<0] <- 0
    survi$conf.type <- method
    return(survi)
}

