#  File R/is.lasttoggle.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
###############################################################################
# is.lasttoggle function tests whether a nw object, with the inherited model 
# info requires duration information. 
###############################################################################
is.lasttoggle <- function(nw,formation=NULL,dissolution=NULL,monitor=NULL,target=NULL){
  
  if(!is.null(formation))
    formation<-ergm.update.formula(formation,nw~., from.new="nw")
  
  if(!is.null(dissolution))  
    dissolution<-ergm.update.formula(dissolution,nw~., from.new="nw")
  
  if(!is.null(monitor)){
    
    unset.offset.call <- function(object){
      if(inherits(object,"call") && object[[1]]=="offset")
        object[[2]]
      else
        object
    }
    
    if(is.character(monitor)){
      monitor <- switch(monitor,
          formation = formation,
          dissolution = dissolution,
          all = append.rhs.formula(~nw, unique(lapply(c(term.list.formula(formation[[3]]),term.list.formula(dissolution[[3]])), unset.offset.call)))
      )
    }
    
    if(!is.null(monitor)) 
      monitor <- ergm.update.formula(monitor,nw~., from.new="nw")
  }
  
  
  if(!is.null(target)){
    if(is.character(targets)){
      targets <- switch(targets,
          formation = formation,
          dissolution = dissolution)}
      
      targets <- ergm.update.formula(targets,nw~., from.new="nw")
    }
  
  
  
    duration.dependent <- if(is.durational(formation) || is.durational(dissolution)|| is.durational(monitor))
        {1} else {0}
    
    duration.dependent
    
  }
