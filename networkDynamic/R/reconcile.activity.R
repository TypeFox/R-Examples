#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2013 the statnet development team
######################################################################

# functions for reconciling edge and vertex activity

reconcile.vertex.activity<-function(net,mode=c("expand.to.edges","match.to.edges","encompass.edges"),edge.active.default=TRUE){
  if (!is.networkDynamic(net)){
    stop("reconcile.vertex.activity can only be applied to networkDynamic objects")
  }
  xn <- substitute(net)

  v<-seq_len(network.size(net))
  mode<-match.arg(mode)
  if(mode=='match.to.edges'){
    # set the activity spells of all vertices to match the activity of their incident edges
    # delete the vertices' activity
    delete.vertex.activity(net,v=v)
    # for each vertex
    for (vert in v) {
      # get the ids of incident edges
      eids<-get.edgeIDs(net,v=vert,neighborhood='combined')
      if (length(eids)>0){
        # get the activity of those edges and  union the edge activiy spells
        # hid the warnings in the case there is no edge activity
        suppressWarnings(activity<-unique(get.edge.activity(net,e=eids,as.spellList=TRUE,active.default=edge.active.default)[,1:2]))
        if (nrow(activity)>0){
          # set vertex activity to edges' activity
          for (i in seq_len(nrow(activity))) {
            activate.vertices(net,v=vert,onset=activity[i,1],terminus=activity[i,2])
          }
        } else {
          # vert's edges have no valid activity, so deactivate
          deactivate.vertices(net,v=vert)
        }
      } else { 
        # vert has no incident edges, so deactivate
        deactivate.vertices(net,v=vert)
      }
    }

  } else if (mode=='expand.to.edges') {
    #vertices activity will be expanded to include the activity periods 
    #of any incident edges (still permits isolated vertices)
    #existing vertex activity is not deleted
    for (vert in v) {
      # get the ids of incident edges
      eids<-get.edgeIDs(net,v=vert,neighborhood='combined')
      if (length(eids)>0){
        # get the activity of those edges and  union the edge activiy spells
        # hide warnings generated if all edgs are inactive
        suppressWarnings(activity<-unique(get.edge.activity(net,e=eids,as.spellList=TRUE,active.default=edge.active.default)[,1:2]))
        if (nrow(activity)>0){
          # set vertex activity to edges' activity
          for (i in seq_len(nrow(activity))) {
            activate.vertices(net,v=vert,onset=activity[i,1],terminus=activity[i,2])
          }
          
        } 
      }
    }
  } else if (mode=='encompass.edges') {
    #vertices activity will be modified so that it has a single spell 
    #beginning with the earliest incident edge activity, and ending with the latest.
    #existing vertex activity is deleted
    delete.vertex.activity(net,v=v)
    for (vert in v) {
      # get the ids of incident edges
      eids<-get.edgeIDs(net,v=vert,neighborhood='combined')
      if (length(eids)>0){
        # get the activity of those edges and  union the edge activiy spells
        activity<-unique(get.edge.activity(net,e=eids,as.spellList=TRUE,active.default=edge.active.default)[,1:2])
        if (nrow(activity)>0){
          # check if the max is a 0-duration spell. if so, include that into the vertex spell
          max.act <- max(activity[,2])
          if (min(activity[,1]) == max.act) {
            # 0 duration activity
          } else {
            # 0 duration spell at the end must be included into vertex spell
            if (any(activity$terminus==max.act & activity$onset==max.act)) max.act <- max.act+1
          }
          
          # set vertex activity to edges' activity
          activate.vertices(net,v=vert,onset=min(activity[,1]),terminus=max.act)
          
        } else {
          deactivate.vertices(net,v=vert)
        } 
      } else {
        deactivate.vertices(net,v=vert)
      }
    }
  } else {
    stop('reconcile vertex activity mode ', mode, ' not supported')
  }
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, net)))
  invisible(net)
}


reconcile.edge.activity<-function(net,mode=c("match.to.vertices","reduce.to.vertices"), active.default=TRUE){
  if (!is.networkDynamic(net)){
    stop("reconcile.edge.activity can only be applied to networkDynamic objects")
  }
  xn <- substitute(net)

  mode<-match.arg(mode)
  v<-seq_len(network.size(net))
  
  if(mode=='reduce.to.vertices'){
    # trim the activity spells of all edges to match the activity of their incident vertices
    
    # for each vertex
    #   get deactivation spells (deacts)
    #   deactivate those spells for incident edges (eids)
    #   reactivate any duration 0 spells (reacts) on those edges (eids)
    
    for (vert in v) {
      eids<-get.edgeIDs(net,v=vert,neighborhood='combined')
      if (length(eids)>0){
        activity <- get.vertex.activity(net,v=vert, as.spellList=TRUE, active.default= active.default)

        if (nrow(activity) > 0) {
          # get the "off" spells
          time.points <- sort(c(activity$onset, activity$terminus))
          
          if (time.points[1] == -Inf) {
            # first active spell is from -Inf
            time.points <- time.points[2:length(time.points)]
          } else {
            # must deactivate from -Inf to first activity point
            time.points <- c(-Inf, time.points)
          }
          if (max(time.points) == Inf) {
            # last active spell is to Inf
            time.points <- time.points[seq_len(length(time.points)-1)]
          } else {
            # must deactivate from last spell to Inf
            time.points <- c(time.points, Inf)
          }
          if (length(time.points) %% 2 != 0) {
            stop(paste("something went wrong with the timeline of vertex", vert))
          }
            
          deacts = matrix(time.points, ncol=2, byrow=T)
    
          # get the 0 duration spells
          activity.0 <- activity$onset[activity$duration == 0]
          
          # need to reactivate the 0 duration vertex spells if they are active for edges
          reacts = lapply(activity.0, function(t) is.active(net, e=eids, at=t))
          
          for (i in seq_len(nrow(deacts))) {
            deactivate.edges(net, e=eids, onset=deacts[i,1], terminus=deacts[i,2])
          }
          
          for (i in seq_along(activity.0)) {
            activate.edges(net, e=eids[reacts[[i]]], at=activity.0[i])
          }
          
        } else {
          # vert has no activity, so deactivate all incident edges
          deactivate.edges(net, e=eids)
        }
      }
      
    }
    
  } else if (mode=='match.to.vertices') {
    # edges will be modified so as to be active whenever all incident vertices are active.
    
    # activate all edges
    # reduce to vertices
    
    activate.edges(net, onset=-Inf, terminus=Inf)
    reconcile.edge.activity(net, mode='reduce.to.vertices', active.default=active.default)
    
  } else {
    stop('reconcile edge activity mode ', mode, ' not supported')
  }
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, net)))
  invisible(net)
}

