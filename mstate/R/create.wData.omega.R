create.wData.omega <-
function(Tstart, Tstop, status, id, stratum, failcode, cens){
  ## Function to create framework for dataset in which individuals are reweighted after they
  ## experience a competing event, the omega-weights in Geskus (2015)
  event.times <- sort(unique(Tstop[status==failcode]))
  n <- length(Tstop)
  Nw <- rep(1,n)
  sel.compet <- (status!=cens)&(status!=failcode)&!is.na(status)
  ## Number of rows in dataset with weights
  Nw[sel.compet] <-  apply(outer(Tstop[sel.compet],event.times,"<"), 1, sum)+1
  data.weight <- data.frame(id=rep(id,Nw), Tstart=NA, Tstop=NA, status=rep(status,Nw),
                            strata=rep(stratum,sum(Nw)))
  data.weight$Tstart <- unlist(lapply(1:n, FUN=function(x,tms,N) {
                                               if (N[x]==1) {
                                                Tstart[x]
                                               } else {
                                                 if (N[x]==2) {
                                                   c(Tstart[x],Tstop[x])
                                                 } else {
                                                   c(Tstart[x], Tstop[x], rev(rev(tms)[2:(N[x]-1)]))
                                                 }
                                               }
                                             }, tms=event.times, N=Nw))
  data.weight$Tstop <- unlist(lapply(1:n, FUN=function(x,tms,N) {
                                              if (N[x]==1) {
                                                Tstop[x]
                                                } else {
                                                  c(Tstop[x], tail(tms,N[x]-1))
                                                }
                                              }, tms=event.times, N=Nw))
  data.weight
}


