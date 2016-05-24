kaplan.meier.location <-
function(fit)
  {
    xloc=fit$time[fit$n.event>0]
    yloc=1-fit$surv[fit$n.event>0]
    nn=length(yloc)
    yloc=(yloc+c(0,yloc[1:(nn-1)]))/2
    return(cbind(xloc,yloc))
  }
