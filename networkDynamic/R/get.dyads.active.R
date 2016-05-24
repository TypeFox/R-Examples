#  File networkDynamic/R/extract.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2012,2013 the statnet development team
######################################################################

# function to return the set of vertex dyads corresponding to active edges in a network
get.dyads.active<-function(nD, onset = NULL, terminus = NULL, length = NULL, at = NULL,  rule = c("any", "all","earliest","latest"), active.default = TRUE){
  
  if(is.hyper(nD)){
    stop("get.dyads.active does not currently support hypergraphic networks")
  }
  rule<-match.arg(rule)
  goodIds<-valid.eids(nD)
  activeEids<-goodIds[is.active(nD,onset=onset,terminus=terminus,length=length,at=at,rule=rule,active.default=active.default,e=goodIds)]
  return(cbind(sapply(nD$mel[activeEids],'[[','outl'),sapply(nD$mel[activeEids],'[[','inl')))
}


#dyadVersion<-function(nw,at){
#  return(get.dyads.active(nw,at=at))
#}

#dfVersion<-function(nw,at){
#  df<-as.data.frame(nw)
#  active<-(df$onset <= at) & (df$terminus >= at)
#  return(cbind(df$tail[active],df$head[active]))  
#}

#times<-microbenchmark(dyadVersion(nw,1),dfVersion(nw,1),times=1000)