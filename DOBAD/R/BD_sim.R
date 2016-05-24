
#
# Forward simulation of process

#note, condCounts 
birth.death.simulant = function(t,X0=1,lambda=1,mu=2,nu=1,
  condCounts=NULL) {
  hx = list()
  hx.times  = c()
  hx.states = c()
  cTime = 0
  cState = X0
  if (is.null(condCounts)){
    while( cTime <= t ) { 
      hx.times  = c(hx.times,cTime)
      hx.states = c(hx.states,cState)
      if( nu==0 && cState == 0 ) {
        break()
      }
      birth = lambda * cState
      death = mu * cState
      rate = nu + birth + death
      cTime = cTime + rexp(n=1,rate=rate)
      event = sample(3,1,prob=c(birth,death,nu)/rate)
      if( event == 1 || event == 3 ) {
        cState = cState + 1
      } else {
        cState = cState - 1
      }
    }
  }
  else {
    while( cTime <= t ) { 
      hx.times  = c(hx.times,cTime)
      hx.states = c(hx.states,cState)
      {#quit if you violate conditional assumptions
        nis <- NijBD(hx.states);
        cNminus <- sum(nis[1,]); #current Nminus
        cNplus <- sum(nis[2,]);
        if (cNminus > condCounts["Nminus"] || cNplus > condCounts["Nplus"]) {
          break();
        }
      }
      if( nu==0 && cState == 0 ) {
        break()
      }
      birth = lambda * cState
      death = mu * cState
      rate = nu + birth + death
      cTime = cTime + rexp(n=1,rate=rate)    
      event = sample(3,1,prob=c(birth,death,nu)/rate)
      if( event == 1 || event == 3 ) {
        cState = cState + 1
      } else {
        cState = cState - 1
      }
    }
  }
  new("BDMC", times=hx.times,states=hx.states,T=t)
}

#
# Calculate summary statistics of a process realization

## count.birth.death.stats = function(simulant) {
##   i = 2
##   cState = simulant$states[1]
##   steps = length(simulant$states)
##   insertions = 0
##   deletions  = 0
##   while( i <= steps ) {
##     if( simulant$states[i] > cState ) {
##       insertions = insertions + 1
##       cState = simulant$states[i]
##     }
##     if( simulant$states[i] < cState ) {
##       deletions = deletions + 1
##       cState = simulant$states[i]
##     }
##     i = i +1
##   }
##   times = c(simulant$times,simulant$T)
##   dwell = sum( diff(times)*simulant$states )
##   return.vector = c(insertions,deletions,cState,dwell)
##   names(return.vector) = c("additions", "removals", "start.state", "particle.timeave")
##   return(return.vector)
## }

#
# Estimate statistics via simulation
#
# n - # of simulants
# lambda - birth rate
# mu - death rate
# nu - immigration rate (=lambda in TFK91)
# X0 - starting state
# marginal - T/F
# end - ending state if marginal=F
#
# Returns:
# c( expected number of births ) - marginal = T
# c( expected number of births ending in state k, probability of ending in state k )
#                  

## simulate.birth.death = function(n=1000,t=1,lambda=1,mu=3, nu=0,end=0,marginal=TRUE,X0=1) {
##   birth = c()
##   death = c()
##   count = c()
##   dwell = c()
##   i = 0
##   total = 0
##   while(i < n) {
##     one = birth.death.simulant(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0)
##     if( marginal ) {
##       stat = count.birth.death.stats(one)
##       birth = c(birth,stat[1])
##       death = c(death,stat[2])
##       count = c(count,stat[1]+stat[2])
##       dwell = c(dwell,stat[4])
##       i = i + 1
##     } else {
##       stat = count.birth.death.stats(one)
##       if( stat[3] == end ) {
##         birth = c(birth,stat[1])
## 	death = c(death,stat[2])
## 	dwell = c(dwell,stat[4])
##         total = total + 1
##       }
##       i = i + 1
##     }
##   }

##   return.vector = NULL
  
##   if( marginal ) {
##     return.vector = c(mean(birth),mean(death),mean(dwell))
##     names(return.vector) = c("add.mean", "rem.mean", "timeave.mean")
##   } else {
##     return.vector = c(total/n, sum(birth)/n, sum(death)/n, sum(dwell)/n,
##       sum(birth)/total, sum(death)/total, sum(dwell)/total)
##     names(return.vector) = c("trans.prob", "add.joint.mean", "rem.joint.mean",
##            "timeave.joint.mean", "add.cond.mean", "rem.cond.mean", "timeave.cond.mean")
##   }
##   return(return.vector)
## }







