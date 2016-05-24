"a.up.low.fun" <-
function(survi, tl, tu)
{
    # Calculates the indices used to derive the critical value
    # determining a Hall-Wellner band. It takes the a survfit object
    # and returns the values belonging to the two timepoints
    # and a matrix with time, the Kaplan-Meier estimator, the sigmas
    # (the sum in Greenwoods formula) and the std. error.
    # t1 should be at least the smallest time, tu the smaller or equal
    # than the highest time.
    
    survi <- survi
    n <- survi$n
    time <- survi$time
    kap.mei <- survi$surv
    indices <- (1:length(time))[survi$n.event>0]
    n.risk <- survi$n.risk
    n.event <- survi$n.event
    a <- n.event/(n.risk*(n.risk-n.event))
    a <- cumsum(a)
    var.st <- kap.mei^2*a    
    std.err <- sqrt(var.st)
    sigma <- var.st/kap.mei^2
    index.low <- max((1:length(time))[(time-tl)<=0])
    index.up <-  max((1:length(time))[(time-tu)<=0])
    sigma.low <- sigma[index.low]
    sigma.up <- sigma[index.up]
    al <- n*sigma.low/(1+n*sigma.low)
    au <- n*sigma.up/(1+n*sigma.up)
    return(list(a.low=al,a.up=au,sigma.mat=cbind(time,kap.mei,sigma,std.err)
           ,start=index.low,end=index.up))
}

