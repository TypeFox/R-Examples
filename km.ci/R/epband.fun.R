"epband.fun" <-
function(survi, tl=NA,tu=NA, method="linear",conf.lev=0.95)
{
    # This function takes a survfit object and modifies it, such that
    # its lower and upper boundaries are now computed using the
    # method by Hall-Wellner.
    # Essentially required are table of critical values,
    # named "critical.value.nair.90", "critical.value.nair.95"
    # "critical.value.nair.99" (see also Appendix C
    # in Klein & Moeschberger p. 451).
    
    data(critical.value.nair.90, critical.value.nair.95, critical.value.nair.99)
    survi <- survi
    tl <- tl
    tu <- tu
    if(max(conf.lev==c(0.90, 0.95, 0.99))!=1)
    {
     stop("confidence level for simultaneous bands must be either 0.90, 0.95 or 0.99")
    }
    # if no tl,tu is given the band covers the whole curve
    if(is.na(tl)&is.na(tu))
    {
        tl <- min(survi$time[survi$n.event>0])
        tu <- max(survi$time[survi$n.event>0 &survi$n.risk>survi$n.event])
    }
    n <- survi$n
    aa <- a.up.low.fun(survi,tl,tu)
    au <- aa$a.up #determines row in table of critical values
    al <- aa$a.low #determines column ...
    dat.mat <- aa$sigma.mat

    # columns used to interpolate
    index.al.left <- floor(al/2*100+1)-1
    index.al.right <- ceiling(al/2*100+1)-1
    if(index.al.left==0)
    {
        index.al.left <- 1
    }
    # rows ...
    index.au.top <- floor((au-0.1)/2*100+1)
    index.au.bottom <- ceiling((au-0.1)/2*100+1)
    if(index.au.bottom==46)
    {
      index.au.bottom <- 45
    }

    #critical values: readingwise 1. topleft,...,3.bottomleft,...
    if(conf.lev==0.90)
    {
    crit1 <- critical.value.nair.90[index.au.top,index.al.left]
    crit2 <- critical.value.nair.90[index.au.top,index.al.right]
    crit3 <- critical.value.nair.90[index.au.bottom,index.al.left]
    crit4 <- critical.value.nair.90[index.au.bottom,index.al.right]
    }
    if(conf.lev==0.95)
    {
    crit1 <- critical.value.nair.95[index.au.top,index.al.left]
    crit2 <- critical.value.nair.95[index.au.top,index.al.right]
    crit3 <- critical.value.nair.95[index.au.bottom,index.al.left]
    crit4 <- critical.value.nair.95[index.au.bottom,index.al.right]
    }
    if(conf.lev==0.99)
    {
    crit1 <- critical.value.nair.99[index.au.top,index.al.left]
    crit2 <- critical.value.nair.99[index.au.top,index.al.right]
    crit3 <- critical.value.nair.99[index.au.bottom,index.al.left]
    crit4 <- critical.value.nair.99[index.au.bottom,index.al.right]
    }
    if(is.na(crit2))# just in case
    {
        crit2 <- (crit1+crit4)/2
    }

    #percentages of interpolation
    vert.perc <- 1-(ceiling(au/2*100)-au/2*100)
    hori.perc <- 1-(ceiling(al/2*100)-al/2*100)

    #interpolations: numbering clockwise
    inter1 <- crit1-(abs(crit1-crit2)*hori.perc)
    inter2 <- crit4-(abs(crit4-crit2)*vert.perc)
    inter3 <- crit3-(abs(crit3-crit4)*hori.perc)
    inter4 <- crit3-(abs(crit3-crit1)*vert.perc)
    interpol <- inter1*(1-vert.perc)+inter4*(1-hori.perc)+inter2*hori.perc+inter3*vert.perc
    crit <- interpol/2

    # First: compute a vector with the deviations
    devia <- abweich.nair.fun(dat.mat,crit)

    # Now, produce a list with the lower and upper boundary
    # dependent of the method.
    if(method=="linear")
    {
        up.low.list <- confi.nair.fun(devia$lin.dev,dat.mat[,2],method)
    }
    if(method=="log")
    {
       up.low.list <- confi.nair.fun(devia$log.dev,dat.mat[,2],method)
    }

    # Finally, modify the survfit object with the new boundaries
    survi$lower <- up.low.list$lower
    survi$upper <- up.low.list$upper
    survi <- modify.surv.fun(survi,aa$start,aa$end,method)

    return(survi)
}

