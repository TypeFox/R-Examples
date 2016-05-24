# functions to return the duration of edge and vertex activity with respect to a query time point


# return the ages of edges activity spell intersecting with at 
edges.age.at <-function(nD,at,e=seq_along(nD$mel),active.default=TRUE){
  if (length(nD$mel) > 0) {
    if ((min(e) < 1) || (max(e) > nD%n% "mnext" - 1)) {
      stop("Illegal edge id in edge.age.at.\n")
    }
  }
  ages<-rep(NA,length(e))
  # find which e are non-deleted edges
  ePresent<-!sapply(nD$mel[e],is.null)
  # restrict to only those present edges that will be active
  ePresent[ePresent]<-is.active(nD,e=e[ePresent],at=at,active.default=active.default)
  # get the spell matrices for those e that are present and active
  activity<-get.edge.activity(nD,e=e[ePresent],active.default=active.default)
  # update age coresponding edges present and active
  ages[ePresent]<-sapply(activity,function(spls){
    index<-spells.hit(c(at,at),spls) # get the index of matching spell
    return(at-spls[index,1]) # return difference between 'at' query and onset of matching spell
  })
  return(ages)
}

# return the ages of  activity spell of edges specified by edge id dyads
dyads.age.at <-function(nD,at,tails,heads,active.default=TRUE,format.out=c('vector','edgelist','matrix')){
  if(is.hyper(nD)){
    stop('dyad.age.at does not support hypergrahic networks')
  }
  if(is.multiplex(nD)){
    stop('dyad.age.at does not support multiplex networks because vertex id dyads may not uniquely specify edges')
  }
  
  format.out<-match.arg(format.out)
  
  # if heads and tails both missing, assume dyads wanted for all existing edges (existing at any time point)
  if(missing(tails) & missing(heads)){
    dyads<-as.matrix.network.edgelist(nD,)
    tails<-dyads[,1]
    heads<-dyads[,2]
  } else {
    if (length(tails)!=length(heads)){
      stop('vectors of vertices for heads and tails of edges must be the same length')
    }
    
  }
  
 
  ages<-rep(NA,length(tails))
  # find the set of edges corresponding to elements in tail and head
  e<-lapply(seq_along(tails),function(i){
    get.edgeIDs.active(nD,v=tails[i],alter=heads[i],at=at,active.default = active.default)
  })
  ePresent<-sapply(e,length)>0
  #if there is at least one active edge..
  if(sum(ePresent)>0){
    activity<-get.edge.activity(nD,e=unlist(e[ePresent]),active.default=active.default)
    # update age coresponding edges present and active
    ages[ePresent]<-sapply(activity,function(spls){
      index<-spells.hit(c(at,at),spls) # get the index of matching spell
      return(at-spls[index,1]) # return difference between 'at' query and onset of matching spell
    })
  }
  if(format.out=='edgelist'){
    out<-cbind(tails,heads,ages)  #todo:  should we force tail head ordering for undirected nets?
  } else if(format.out=='matrix'){
    out<-matrix(NA,nrow=network.size(nD),ncol=network.size(nD))
    
    #I'm not sure why this is the inexing form but,  out[tails,heads] doesn't work
    out[tails+(heads-1)*network.size(nD)]<-ages 
    # if it is undirected, also set the other dimension
    if(!is.directed(nD)){
      out[heads+(tails-1)*network.size(nD)]<-ages 
    }
    
  } else { # return as a vector
    out<-ages
  }
  
  return(out)
  
  
}



vertices.age.at <-function(nD,at,v=seq_len(network.size(nD)),active.default=TRUE){
  ages<-rep(NA,length(v))
  vActive<-is.active(nD,v=v,at=at,active.default=active.default)
  activity<-get.vertex.activity(nD,v=v[vActive],active.default=active.default)
  # update age coresponding to vertices  active
  ages[vActive]<-sapply(activity,function(spls){
    index<-spells.hit(c(at,at),spls) # get the index of matching spell
    return(at-spls[index,1]) # return difference between 'at' query and onset of matching spell
  })
  return(ages)
}