# functions for determining when network elements take specific values


# internal function to determine the first / earliest time point at which a tea attribute matches using the specified operation
whenMatchEarliest<-function(tea,value,match.op='==',no.match=Inf){
  sapply(tea,function(val){
    if (length(val)==1 && is.na(val) | is.null(val)){
      return(NA)
    } else {
      #value list
      comparisons<-sapply(val[[1]],match.op,value)
      if(!is.logical(comparisons)){
        stop("the comparison function provided by the 'match.op' argument must provide Logical results for every attribute value")
      }
      matches<-which(comparisons)
      if(length(matches)<1){
        return(no.match)
      } 
      return(val[[2]][min(matches),1])
    }
  })
}

# internal function to determine the last / latest time point at which a tea attribute matches using the specified operation
whenMatchLatest<-function(tea,value,match.op='==',no.match=Inf){
  sapply(tea,function(val){
    if (length(val)==1 && is.na(val) | is.null(val)){
      return(NA)
    } else {
      #value list
      comparisons<-sapply(val[[1]],match.op,value)
      if(!is.logical(comparisons)){
        stop("the comparison function provided by the 'match.op' argument must provide Logical results for every attribute value")
      }
      matches<-which(comparisons)
      if(length(matches)<1){
        return(no.match)
      } 
      return(val[[2]][max(matches),2])
    }
  })
}

# function to determine when a TEA vertex attribute takes on a specific value

when.vertex.attrs.match<-function(nd,attrname,value,match.op='==',rule='earliest',no.match=Inf,v=seq_len(network.size(nd))){
  if(missing(nd) || !is.networkDynamic(nd)){
    stop('when.vertex.attrs.match requires its first argument to be a networkDynamic object')
  }
  if(missing(value)){
    stop("when.vertex.attrs.match requires providing a 'value' argument for comparison")
  }
  if(missing(attrname)){
    stop("when.vertex.attrs.match requires providing an 'attrname' argument giving the name of the vertex attribute to be compared")
  }
  
  # search for the .active version of the name
  searchAttr<-paste(attrname,"active",sep='.')
  tea<-get.vertex.attribute(nd,searchAttr,unlist=FALSE)[v]
  # use different method depending on rule
  if (rule=='earliest'){
    return(whenMatchEarliest(tea=tea,value=value,match.op=match.op,no.match=no.match))
  } else if (rule=='latest'){
    return(whenMatchLatest(tea=tea,value=value,match.op=match.op,no.match=no.match))
  } else {
    stop("no matching methods implemented for rule '",rule,"'")
  }

}


# function to determine when a TEA edge attribute takes on a specific value

when.edge.attrs.match<-function(nd,attrname,value,match.op='==',rule='earliest',no.match=Inf,e=seq_along(nd$mel)){
  if(missing(nd) || !is.networkDynamic(nd)){
    stop('when.vedge.attrs.match requires its first argument to be a networkDynamic object')
  }
  if(missing(value)){
    stop("when.edge.attrs.match requires providing a 'value' argument for comparison")
  }
  if(missing(attrname)){
    stop("when.edge.attrs.match requires providing an 'attrname' argument giving the name of the edge attribute to be compared")
  }
  
  # search for the .active version of the name
  searchAttr<-paste(attrname,"active",sep='.')
  tea<-get.edge.attribute(nd,searchAttr,unlist=FALSE)[e]
  # use different method depending on rule
  if (rule=='earliest'){
    return(whenMatchEarliest(tea=tea,value=value,match.op=match.op,no.match=no.match))
  } else if (rule=='latest'){
    return(whenMatchLatest(tea=tea,value=value,match.op=match.op,no.match=no.match))
  } else {
    stop("no matching methods implemented for rule '",rule,"'")
  }
}


# function to find the next time at which an edge involving the specified subset of vertices is activated or deactivated
when.next.edge.change<-function(nd,at, v=seq_len(network.size(nd)),neighborhood = c("out", "in", "combined")){
  # find ids of the set of edges we need to check
  eids<-unique(unlist(sapply(v,function(v){
    get.edgeIDs(nd,v=v,neighborhood=neighborhood) 
  })))
  # find which ones are currently active
  nextTimes<-sapply(eids,function(e){
    # if no activitey return inf
    spls<-nd$mel[[e]]$atl$active
    if (is.null(spls)){
      return(Inf)
    } 
    splIndex<-spells.hit(needle=c(at,Inf),haystack=spls)
    if (splIndex<0){
      return(Inf)
    }
    # if the spell is active, return its terminus
    if(spls[splIndex,1]<=at){
      return(spls[splIndex,2])
    } else {
      # otherwise return its onset
      return(spls[splIndex,1])
    }
  })
  return(min(nextTimes))
}