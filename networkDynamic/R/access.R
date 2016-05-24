#  File networkDynamic/R/access.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
# This file contains various routines for accessing network class objects with
# dynamic extensions.
#
# Contents:
#
#   activate.edges
#   activate.vertices
#   add.edges.active
#   add.vertices.active
#   deactivate.edges
#   deactivate.vertices
#   delete.edge.activity
#   delete.vertex.activity
#   get.edgeIDs.active
#   get.edges.active
#   get.neighborhood.active
#   get.change.times
#   is.active
#   is.adjacent.active
#   network.dyadcount.active
#   network.edgecount.active
#   network.naedgecount.active
#   network.size.active
#   insert.spell
#   delete.spell
#
######################################################################

#Function to activate the selected edges at the appropriate time.  If already
#active, activation has no effect; otherwise, it inserts an onset time at
#the appropriate mark.  Edges without an "active" attribute are given one.
activate.edges <- function(x, onset=NULL, terminus=NULL, length=NULL, at=NULL,
                           e=seq_along(x$mel)){
  xn <- substitute(x)   # needed for proper assignment in calling environment

  # checks for proper inputs, translations into onset and terminus
  if(!is.network(x)) 
    stop("activate.edges requires an argument of class network.\n")
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Activation times must be a numeric vector in activate.edges.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in activate.edges.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in activate.edges.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in activate.edges.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.vector(e) || !is.numeric(e))
    stop("Edge ID's, e, must be a numeric vector in activate.edges.\n")
  if(length(x$mel)>0) {
    if((min(e,Inf) < 1) || (max(e,-Inf) > x%n%"mnext"-1)) 
      stop("Illegal edge in activate.edges.\n")
 
    # preliminaries
    e <- e[!sapply(x$mel[e], is.null)]  #Filter out non-edges
    if(!is.null(at)) {
      onset <- terminus <- rep(at, length=length(e))
    } else if (!is.null(onset)) {
      onset <- rep(onset, length=length(e))
      if(!is.null(terminus))
        terminus <- rep(terminus, length=length(e))
      else if (!is.null(length))
        terminus <- onset + rep(length, length=length(e))
    } else {
      if (is.null(terminus)) {
        onset <- rep(-Inf, length=length(e))
        terminus <- rep(Inf, length=length(e))
      } else {
        terminus <- rep(terminus, length=length(e))
        onset <- terminus - rep(length, length=length(e))
      }
    }
    if(any(onset>terminus))
      stop("Onset times must precede terminus times in activate.edges.\n")
    
    x <- .Call(ActivateEdges_R, x, onset, terminus, e, FALSE)
  }
  
  set.nD.class(x)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)  
}


#Function to activate the selected vertices at the appropriate time.  If already
#active, activation has no effect; otherwise, it inserts an onset time at
#the appropriate mark.  Vertices without an "active" attribute are given one.
activate.vertices <- function(x, onset=NULL, terminus=NULL, length=NULL, at=NULL,
                              v=seq_len(network.size(x))) {
  xn <- substitute(x)   # needed for proper assignment in calling environment
  
  # checks for proper inputs
  if(!is.network(x)) 
    stop("activate.vertices requires an argument of class network.\n")
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Activation times must be a numeric vector in activate.vertices.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in activate.vertices.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in activate.vertices.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in activate.vertices.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.vector(v) || !is.numeric(v))
    stop("Vertex ID's, v, must be a numeric vector in activate.vertices.\n")
  if((min(v,Inf) < 1) || (max(v,-Inf) > network.size(x))) 
    stop("Illegal vertex in activate.vertices.\n")

  # preliminaries
  v <- v[!sapply(x$val[v], is.null)]  #Filter out non-vertices
  if(!is.null(at)) {
    onset <- terminus <- rep(at, length=length(v))
  } else if (!is.null(onset)) {
    onset <- rep(onset, length=length(v))
    if(!is.null(terminus))
      terminus <- rep(terminus, length=length(v))
    else if (!is.null(length))
      terminus <- onset + rep(length, length=length(v))
  } else {
    if (is.null(terminus)) {
      onset <- rep(-Inf, length=length(v))
      terminus <- rep(Inf, length=length(v))
    } else {
      terminus <- rep(terminus, length=length(v))
      onset <- terminus - rep(length, length=length(v))
    }
  }
  if(any(onset>terminus))
    stop("Onset times must precede terminus times in activate.vertices.\n")
  
  # choosing to ignore activation requests of (Inf,Inf) or (-Inf, -Inf)
  ignore <- (onset==Inf) | (terminus==-Inf)
  if(any(ignore)){
    onset<-onset[!ignore]; terminus<-terminus[!ignore]; v<-v[!ignore]
  }
  if(length(v) > 0) {

    # get current active matrices and insert spells
    uniqueV<-unique(v)
    active <- lapply(x$val[uniqueV], "[[", "active")
    infMat<-matrix(c(-Inf,Inf),1,2)
    for(i in 1:length(v)){
      index<-which(uniqueV==v[i])
      if(!(identical(active[[index]], infMat)))
        active[[index]] <- insert.spell(active[[index]], onset[i], terminus[i])
    }
    set.vertex.attribute(x, "active", active, uniqueV)
  }

  set.nD.class(x)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)  
}


#Function to deactivate the selected edges at the appropriate time.  If already
#inactive, activation has no effect; otherwise, it inserts a termination time at
#the appropriate mark.  Edges without an "active" attribute are given one.
deactivate.edges<-function(x, onset=NULL, terminus=NULL, length=NULL, at=NULL,
                           e=seq_along(x$mel)){
  xn <- substitute(x)   # needed for proper assignment in calling environment

  # checks for proper inputs
  if(!is.network(x)) 
    stop("deactivate.edges requires an argument of class network.\n")
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Deactivation times must be a numeric vector in deactivate.edges.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in deactivate.edges.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in deactivate.edges.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in deactivate.edges.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.vector(e) || !is.numeric(e))
    stop("Edge ID's, e, must be a numeric vector in deactivate.edges.\n")
  if(length(x$mel) > 0) {
    if((min(e,Inf) < 1) || (max(e,-Inf) > x%n%"mnext"-1)) 
      stop("Illegal edge in deactivate.edges.\n")

    # preliminaries
    e <- e[!sapply(x$mel[e], is.null)]  #Filter out non-edges
    if(length(e)==0) return(invisible(set.nD.class(x)))
    if(!is.null(at)) {
      onset <- terminus <- rep(at, length=length(e))
    } else if (!is.null(onset)) {
      onset <- rep(onset, length=length(e))
      if(!is.null(terminus))
        terminus <- rep(terminus, length=length(e))
      else if (!is.null(length))
        terminus <- onset + rep(length, length=length(e))
    } else {
      if (is.null(terminus)) {
        onset <- terminus <- rep(Inf, length=length(e))
      }else {
        terminus <- rep(terminus, length=length(e))
        onset <- terminus - rep(length, length=length(e))
      }
    }
    if(any(onset>terminus))
      stop("Onset times must precede terminus times in deactivate.edges.\n")

  
    #Get existing activity attributes and update as needed
    active<-lapply(lapply(x$mel[e],"[[","atl"),"[[","active")
    for(i in seq_along(active)){
      if(is.infinite(onset[i]) && is.infinite(terminus[i])){
        active[[i]]<-matrix(c(Inf, Inf),1,2)
      }else if(is.null(active[[i]])){
        if(is.infinite(onset[i]))
          active[[i]]<-matrix(c(terminus[i],Inf),1,2)
        else if (is.infinite(terminus[i]))
          active[[i]]<-matrix(c(-Inf,onset[i]),1,2)
        else
          active[[i]]<-matrix(c(-Inf,terminus[i],onset[i],Inf),2,2)
      }else if(!all(active[[i]]==Inf) && !all(active[[i]]==-Inf)){
        active[[i]] <- delete.spell(active[[i]], onset[i], terminus[i])
      }
    }
    set.edge.attribute(x=x,attrname="active",value=active,e=e)
  }
  
  set.nD.class(x)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)  
}


#Function to deactivate the selected vertices at the appropriate time.  If 
#already inactive, activation has no effect; otherwise, it inserts a termination
#time at the appropriate mark.  Vertices without an "active" attribute are given
#one.
deactivate.vertices<-function(x, onset=NULL, terminus=NULL, length=NULL, at=NULL,
                               v=seq_len(network.size(x)), deactivate.edges=FALSE){
  xn <- substitute(x)   # needed for proper assignment in calling environment

  # checks for proper inputs
  if(!is.network(x)) 
    stop("deactivate.vertices requires an argument of class network.\n")
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Deactivation times must be a numeric vector in deactivate.vertices.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in deactivate.vertices.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in deactivate.vertices.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in deactivate.vertices.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.vector(v) || !is.numeric(v))
    stop("Vertices, v, must be a numeric vector in deactivate.vertices.\n")
  
  if((min(v,Inf) < 1) || (max(v,-Inf) > network.size(x))) 
    stop("Illegal vertex in deactivate.vertices.\n")
  
  # preliminaries
  v <- v[!sapply(x$val[v], is.null)]  #Filter out non-vertices
  if(length(v) > 0) {
    if(!is.null(at)) {
      onset <- terminus <- rep(at, length=length(v))
    } else if (!is.null(onset)) {
      onset <- rep(onset, length=length(v))
      if(!is.null(terminus))
        terminus <- rep(terminus, length=length(v))
      else if (!is.null(length))
        terminus <- onset + rep(length, length=length(v))
    } else {
      if (is.null(terminus)) {
        onset <- terminus <- rep(Inf, length=length(v))
      }else {
        terminus <- rep(terminus, length=length(v))
        onset <- terminus - rep(length, length=length(v))
      }
    }
    if(any(onset>terminus))
      stop("Onset times must precede terminus times in deactivate.vertices.\n")

    #Get existing activity attributes and update as needed
    active<-lapply(x$val[v],"[[","active")
    for(i in seq_along(active)){
      if(is.infinite(onset[i]) && is.infinite(terminus[i])){
        active[[i]]<-matrix(c(Inf, Inf),1,2)
      }else if(is.null(active[[i]])){
        if(is.infinite(onset[i]))
          active[[i]]<-matrix(c(terminus[i],Inf),1,2)
        else if (is.infinite(terminus[i]))
          active[[i]]<-matrix(c(-Inf,onset[i]),1,2)
        else
          active[[i]]<-matrix(c(-Inf,terminus[i],onset[i],Inf),2,2)
      }else if(!all(active[[i]]==Inf) && !all(active[[i]]==-Inf)){
        active[[i]] <- delete.spell(active[[i]], onset[i], terminus[i])
      }
    }
    set.vertex.attribute(x=x,attrname="active",value=active,v=v)
  }
  
  # deactivate the associated edges, if user wants
  if (deactivate.edges) {
    e = NULL
    for (vert in v) {
      e = c(e, get.edgeIDs.active(x, v=vert, onset=onset, terminus=terminus, 
                                  length=length, at=at, neighborhood="combined"))
    }
    if (length(e) > 0) {
      deactivate.edges(x, onset=onset, terminus=terminus, length=length, at=at, e=unique(e))
    }
  }
  
  set.nD.class(x)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)  
}


# adds new edges, active at the given time
add.edges.active <- function(x, tail, head, names.eval=NULL, vals.eval=NULL, onset=NULL, terminus=NULL, length=NULL, at=NULL, ...) {
  xn <- substitute(x)   # needed for proper assignment in calling environment

  if(!is.network(x))
    stop("add.edges.active requires an argument of class network.\n")
  if(!is.numeric(tail) || !is.numeric(head))
    stop("The vertex ID's given in 'tail' and 'head' must be a numeric vector in add.edges.active.\n")
  if(min(tail) < 1 || max(tail) > network.size(x))
    stop("Illegal vertex in 'tail' vector in add.edges.active .\n")
  if(min(head) < 1 || max(head) > network.size(x))
    stop("Illegal vertex in 'head' vector in add.edges.active .\n")

  n = max(length(tail), length(head))
  if(length(tail) != length(head)) {
    tail = rep(tail, length=n)
    head = rep(head, length=n)
  }

  add.edges(x, tail, head,names.eval,vals.eval)
  activate.edges(x, onset, terminus, length, at, e=seq(x%n%"mnext"-n, x%n%"mnext"-1))
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)  
}


# adds new vertices, active at the given time
add.vertices.active <- function(x, nv, vattr=NULL, last.mode=TRUE, onset=NULL, terminus=NULL, length=NULL, at=NULL,...) {
 
  if(!is.network(x))
    stop("add.vertices.active requires an argument of class network.\n")
  if(!is.numeric(nv))
    stop("The number of vertices given in 'nv' must be numeric in add.verices.active.\n")
  xn <- substitute(x)   # needed for proper assignment in calling environment

  if (nv>0){
    add.vertices(x, nv,vattr,last.mode)
    activate.vertices(x, onset, terminus, length, at, v=seq(x%n%"n"-nv+1, x%n%"n"))
  } else {
    if(!is.networkDynamic(x)){
      x<-set.nD.class(x)
    }
  }
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)  
}

# ------------- get.change.times ---------
# pulls out all of the times at which acitvity changes
# TODO: may be problems with the 'null' (Inf,Inf) spell
get.change.times <- function (x, vertex.activity=TRUE,edge.activity=TRUE, ignore.inf=TRUE,vertex.attribute.activity=TRUE,edge.attribute.activity=TRUE,network.attribute.activity=TRUE) {
  if(!is.network(x))
    stop("get.change.times requires an argument of class network.\n")
  if(!is.logical(vertex.activity) | !is.logical(edge.activity))
    stop("get.change.times requires that vertex.activity and edge.activity be logicals.\n")
  times <- numeric(0)
  if(edge.activity & network.edgecount(x)>0){
    spls<-get.edge.attribute(x$mel, "active",unlist=FALSE)
    spls<-.removeNullSpells(spls)
    times <- c(times,unlist(spls) )
  }
  if(vertex.activity & network.size(x)>0){
    if("active"%in%list.vertex.attributes(x)){
      spls<-get.vertex.attribute(x, "active",unlist=FALSE)
      spls<-.removeNullSpells(spls)
      times <- c(times, unlist(spls))
    }
  }
  if(vertex.attribute.activity & network.size(x)>0){
    attrs<-list.vertex.attributes.active(x,onset=-Inf,terminus=Inf,dynamic.only=TRUE)
    for(attr in attrs){
      vals <- get.vertex.attribute.active(x,sub('.active','',attr),onset=-Inf,terminus=Inf,return.tea=TRUE)
      vals <- vals[!is.na(vals)]
      times<-c(times, unique(unlist(sapply(vals,'[[',2,simplify=FALSE))))
    }
  }
  if(edge.attribute.activity & network.edgecount(x)>0){
    attrs<-list.edge.attributes.active(x,onset=-Inf,terminus=Inf,dynamic.only=TRUE)
    for(attr in attrs){
      vals<-get.edge.attribute.active(x,sub('.active','',attr),onset=-Inf,terminus=Inf,return.tea=TRUE)
      vals<-vals[!is.na(vals)]       
      times<-c(times, unique(unlist(sapply(vals,'[[',2,simplify=FALSE))))
    }
  }
  if(network.attribute.activity){
    attrs<-list.network.attributes.active(x,onset=-Inf,terminus=Inf,dynamic.only=TRUE)
    for(attr in attrs){
      times<-c(times, unique(as.vector(get.network.attribute.active(x,sub('.active','',attr),onset=-Inf,terminus=Inf,return.tea=TRUE)[[2]])))
    }
  }
  
  if(ignore.inf){
    times <- sort(unique(times[!is.infinite(times)]))
  } else {
    times <- sort(unique(times))
  }
  return(times)
}


#Variant of get.edgeIDs with dynamic query support
get.edgeIDs.active<-function(x,v,onset=NULL,terminus=NULL,length=NULL, at=NULL,
                             alter=NULL,neighborhood=c("out", "in", "combined"),
                             rule=c("any","all","earliest","latest"),na.omit=TRUE,active.default=TRUE){
  
  if(missing(v)){
    stop("'v' parameter must be specified with a vertex id to indicate which vertex to search for incident edges")
  }
  rule<-match.arg(rule)
  # get IDs and filter by activity
  eid<-get.edgeIDs(x=x,v=v,alter=alter,neighborhood=neighborhood, na.omit=na.omit)
  if(length(eid)==0)
    return(integer(0))
  active = is.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
    e=eid,v=NULL,rule=rule, active.default=active.default)
  if(!any(active))
    return(integer(0))
  eid[active]
}


#Variant of get.edges with dynamic query support.  (Note: not safe in the long
#run...)
get.edges.active<-function(x,v,onset=NULL,terminus=NULL,length=NULL, at=NULL,
                           alter=NULL,neighborhood=c("out", "in", "combined"),
                           rule=c("any","all","earliest","latest"),na.omit=TRUE,active.default=TRUE){
  if(missing(v)){
    stop("'v' parameter must be specified with vertex id to indicate which vertex to search for incident edges")
  }
  rule<-match.arg(rule)
  # get IDs and filter by activity
  eid<-get.edgeIDs(x=x,v=v,alter=alter,neighborhood=neighborhood, na.omit=na.omit)
  if(length(eid)==0)
    return(list())
  active = is.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
    e=eid,v=NULL,rule=rule, active.default=active.default)
  if(!any(active))
    return(list())
  x$mel[eid][active]
}


#Variant of get.neighborhood with dynamic query support.  Slow, most likely.
get.neighborhood.active<-function(x,v,onset=NULL,terminus=NULL,length=NULL, at=NULL,
                                  type=c("out", "in", "combined"),rule=c("any","all","earliest","latest"),na.omit=TRUE,active.default=TRUE){
  rule<-match.arg(rule)
  # get active edges and assemble neighborhood in questions
  if(!is.directed(x)){
    el<-get.edges.active(x=x,v=v,onset=onset,terminus=terminus,length=length, at=at,
                         alter=NULL, neighborhood="out", rule=rule,na.omit=na.omit,
                         active.default=active.default)
    if(length(el)>0){
      neigh<-sort(unique(c(sapply(el,"[[","inl"),sapply(el,"[[","outl"))))
      #Loop check
      if(!any(sapply(el,function(z){(v%in%z[["inl"]])&&(v%in%z[["outl"]])})))
        neigh<-neigh[neigh!=v]
    }else
      neigh<-integer(0)
  }else{   
    if(match.arg(type)=="out"){   # directed out neighboorhood
      el<-get.edges.active(x=x,v=v,onset=onset,terminus=terminus,length=length, at=at,
                           alter=NULL, neighborhood="out",rule=rule,na.omit=na.omit,
                           active.default=active.default)
      if(length(el)>0)
        neigh<-sort(unique(sapply(el,"[[","inl")))
      else
        neigh<-integer(0)
    }else if(match.arg(type)=="in"){   # directed in neighboorhood
      el<-get.edges.active(x=x,v=v,onset=onset,terminus=terminus,length=length, at=at,
                           alter=NULL, neighborhood="in",rule=rule,na.omit=na.omit,
                           active.default=active.default)
      if(length(el)>0){
        neigh<-sort(unique(sapply(el,"[[","outl")))
      }else
        neigh<-integer(0)
    }else{                            # directed in/out neighboorhood
      out.el<-get.edges.active(x=x,v=v,onset=onset,terminus=terminus,length=length,at=at,
                               alter=NULL, neighborhood="out",rule=rule,na.omit=na.omit,
                               active.default=active.default)
      if(length(out.el)>0)
        neigh<-sort(unique(sapply(out.el,"[[","inl")))
      else
        neigh<-integer(0)
      in.el<-get.edges.active(x=x,v=v,onset=onset,terminus=terminus,length=length,at=at,
                              alter=NULL, neighborhood="in",rule=rule,na.omit=na.omit,
                              active.default=active.default)
      if(length(in.el)>0)
        neigh<-sort(unique(c(neigh,sapply(in.el,"[[","outl"))))
    }
  }
  neigh
}

# wrapper functions to return activity matrices of edges and vertices
# THESE WERE NOT BEING USED, see version in utilities.R


#Function to assess activity of edges (e) or vertices (v) at a given point
#or in a given interval.  If an interval is specified, then rule=="any" 
#returns TRUE for elements active at any time in the interval.  The rule=="all"
#setting returns TRUE for elements active during the entire interval.  Unless
#given either e or v, the function returns NA.
#
#Note that there are a lot of complications here surrounding Inf values.  If
#an activity spell starts at time Inf, it can never match anything (including
#query onsets of Inf).  If an activity spell starts at finite time and ends
#at Inf, however, it _does_ match an onset/terminus of Inf.  By turns, a 
#spell which begins at time -Inf should match -Inf onset times.  All this is
#very annoying, and makes me wish that I'd just outlawed infinity.  But that's
#how things are.
is.active<-function(x,onset=NULL,terminus=NULL,length=NULL, at=NULL, e=NULL,v=NULL, rule=c("any","all","earliest","latest"),active.default=TRUE){
  # checks for proper inputs
  if(!is.network(x)) 
    stop("is.active requires an argument of class network.\n")
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Singular time points given by the 'at' argument must be a numeric vector in is.active.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in is.active.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector is.active.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in is.active.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(is.null(terminus) || is.null(length))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(length(e)*length(v)>0)
    stop("Either edges or vertices must be specified (not both) in is.active.\n")
  if(!is.null(v) & length(v)>0){
    if(!is.vector(v) || !is.numeric(v))
      stop("Vertex ID's, v, must be a numeric vector in is.active.\n")
    if((min(v) < 1) || (max(v) > network.size(x))) 
      stop("Vertex ID's, v, must be in the range from 1 to the size of the network in is.active.\n")
  }
  if(!is.null(e)){
    if(!is.vector(e) || !is.numeric(e))
      stop("Edge ID's, e, must be a numeric vector in is.active.\n")
    if((min(e,Inf) < 1) || (max(e,-Inf) > x%n%"mnext"-1)) 
      stop("Edge ID's in is.active e argument must be in the range from 1 to the number of edges in the network.\n")
  }

  # vertices or edges?
  if(length(e)){
    origelen<-length(e)
    e <- e[!sapply(x$mel[e], is.null)]  # filter out non-edges
    # if e were omitted due to null edges, give warning
    if (length(e)< origelen){
      warning("Some edge IDs in the e argument correspond to deleted edges and will be ignored. Indices of values returned will not correspond to elements of e.")
    }
  }
  if(length(v))
    v <- v[!sapply(x$val[v], is.null)]  # filter out non-vertices TODO: can this happen?
  if(length(e)+length(v)==0)
    return(logical(0))
  if(length(e)){
    active<-lapply(lapply(x$mel[e],"[[","atl"),"[[","active")
    ev <- e  
  } else { 
    active<-lapply(x$val[v],"[[","active")
    ev <- v
  }

  # preliminaries
  rule<-match.arg(rule)
  if(!is.null(at)) {
    onset <- terminus <- rep(at, length=length(ev))
  } else if (!is.null(onset)) {
    onset <- rep(onset, length=length(ev))
    if(!is.null(terminus))
      terminus <- rep(terminus, length=length(ev))
    else if (!is.null(length))
      terminus <- onset + rep(length, length=length(ev))
  } else {
    terminus <- rep(terminus, length=length(ev))
    onset <- terminus - rep(length, length=length(ev))
  }
  if(any(onset>terminus))
    stop("Onset times must precede terminus times in is.active.\n")

#  return(.Call('IsActiveInVector', onset, terminus, active, (match.arg(rule) == 'all'), active.default, get("debug.output", envir=.GlobalEnv)))
  return(.Call(IsActiveInVector_R, onset, terminus, active, (match.arg(rule) == 'all'), active.default, FALSE))
}


#Variant of is.adjacent for networks with dynamic extensions.  Slow, but will
#get the job done.
is.adjacent.active<-function(x,vi,vj,onset=NULL,terminus=NULL,length=NULL, at=NULL,
                             rule=c("any","all","earliest","latest"),na.omit=FALSE,active.default=TRUE){
  rule<-match.arg(rule)
  #Initially, get edge IDs from vi to vj
  eid<-get.edgeIDs(x=x,v=vi,alter=vj,neighborhood="out",na.omit=na.omit)
  #Return TRUE iff any active edges exist
  if(length(eid)==0)
    FALSE
  else
    any(is.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
                  e=eid,v=NULL,rule=rule, active.default=active.default))
}


#Variant network.dyadcount which uses only active vertices.
network.dyadcount.active<-function (x, onset=NULL, terminus=NULL, length=NULL, at=NULL,
                                    rule=c("any","all","earliest","latest"), na.omit = TRUE, active.default=TRUE) {
  if (!is.network(x)) 
    stop("network.dyadcount.active requires an argument of class network.")
  rule<-match.arg(rule)
  if (is.bipartite(x)) {
    bip = x%n%"bipartite"
    nactor <- ifelse(bip >= 0,sum(is.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
                                           v=seq_len(bip),rule=rule, active.default=active.default)),
                     0)
    nevent <- ifelse(x%n%"n">0 && bip<x%n%"n",sum(is.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
                                                          v=(bip+1):(x%n%"n"),rule=rule, active.default=active.default)),
                     0)
    if (is.directed(x))
      dyads <- nactor * nevent * 2
    else
      dyads <- nactor * nevent
  } else {
    nodes <- ifelse(x%n%"n">0, network.size.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
                                                   rule=rule, active.default=active.default),
                    0)
    if (is.directed(x)) 
      dyads <- nodes * (nodes - 1)
    else 
      dyads <- nodes * (nodes - 1)/2
  }

  if (na.omit && x%n%"mnext" > 1 && dyads > 0) {
    # note that I've removed a code block that would replace the block below.
    # it handles the count of missing edges through the design attribute, rather
    # than through the network.naedgecount.active function. The code chunk
    # can be found in the v0.1 tag.  --alc

    # second note, you cannot just count up the number of na edges and subtract
    # this from 'dyads', since 'dyads' is counted over the smaller active network,
    # and 'na.edgecount' counts active missing edges over the full network.
    # you can have edges that are active, whose head/tail nodes are not, and this
    # leads to incorrect counts.
    # given that we're supposed to upload to cran tomorrow, I'm using the quick
    # fix:  this will be slow.  We (I) should fix this later.
    
    xextracted = network.extract(x=x,onset=onset,terminus=terminus,length=length,at=at,
      rule=rule,active.default=active.default)
    na.edgecount = network.naedgecount.active(x=xextracted,onset=onset,terminus=terminus,
      length=length,at=at,rule=rule, active.default=active.default)
    dyads <- dyads - na.edgecount
  }
  dyads
}


#Variant network.edgecount which counts only active edges.  Not long-run safe.
network.edgecount.active<-function (x, onset=NULL, terminus=NULL, length=NULL, at=NULL,
                                    rule=c("any","all","earliest","latest"), na.omit = TRUE, active.default=TRUE){
  rule<-match.arg(rule)
  if(x%n%"mnext">1){
    act<-is.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
                 e=valid.eids(x), v=NULL,rule=rule, active.default=active.default)
    if(na.omit)
      sum(act*(1-(x%e%"na")))
    else
      sum(act)
  } else {
    0
  }
}


#Variant network.naedgecount which counts only active edges.  Not safe.
network.naedgecount.active<-function (x, onset=NULL, terminus=NULL, length=NULL, at=NULL,
                                      rule=c("any","all","earliest","latest"), active.default=TRUE){
  if(x%n%"mnext">1) {
    act<-is.active(x=x,onset=onset,terminus=terminus,length=length, at=at,
                   e=valid.eids(x),v=NULL,rule=rule, active.default=active.default)
    sum(act*(x%e%"na"))
  }else{
    0
  }
}

#Network size which counts only active vertices - don't use for other purposes!
network.size.active<-function(x,onset=NULL,terminus=NULL,length=NULL, at=NULL,
                              rule=c("any","all","earliest","latest"),active.default=TRUE){
  rule<-match.arg(rule)
  sum(is.active(x=x,onset=onset,terminus=terminus,length=length, at=at,
                e=NULL,v=seq_len(network.size(x)), rule=rule,active.default=active.default))
}


#--------------------------------------------------------------
# this function removes the activity matrices for a given
# set of edges.
#
# @param
#    x: a networkDynamic or network object
#    e: the edges whose spell matrices are to be deleted;
#       default=all
#
# @return:
#    the networkDynamic object without the spell matrices of 'e'
#------------------------------------------------------------------
delete.edge.activity <- function(x, e=seq_along(x$mel)) {
  xn <- substitute(x)   # needed for proper assignment in calling environment
  
  if(!is.network(x))
    stop("The remove.activity function requires that x be a network object.\n")
  if(!is.vector(e) || !is.numeric(e))
    stop("Edge ID's, e, must be a numeric vector in remove.activity.\n")
  if((min(e,Inf) < 1) || (max(e,-Inf) > x%n%"mnext"-1)) 
    stop("Illegal edge in remove.activity argument e.\n")

  # if deleting all edges, can use network's delete.edge.attribute
  # function, otherwise need to manually remove activity matrices.
  if (length(e) == length(x$mel)) {
    delete.edge.attribute(x, "active")
  } else {
    leave.active = setdiff(seq_along(x$mel), e)
    leave.active = leave.active[!sapply(x$mel[leave.active], is.null)]  # filter out non-edges
    left.activity = lapply(lapply(x$mel[leave.active], "[[", "atl"), "[[", "active")
    delete.edge.attribute(x, "active")
    set.edge.attribute(x, "active", left.activity, leave.active)
  }
  
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)  
}


#--------------------------------------------------------------
# this function removes the activity matrices for a given
# set of vertices.
#
# @param
#    x: a networkDynamic or network object
#    v: the vertices whose spell matrices are to be deleted;
#       default=all
#
# @return:
#    the networkDynamic object without the spell matrices of 'v'
#------------------------------------------------------------------
delete.vertex.activity <- function(x, v=seq_len(network.size(x))) {
  xn <- substitute(x)   # needed for proper assignment in calling environment
  
  if(!is.network(x))
    stop("The remove.activity function requires that x be a network object.\n")
  if(!is.vector(v) || !is.numeric(v))
    stop("Vertex ID's, v, must be a numeric vector in remove.activity.\n")
  if((min(v,Inf) < 1) || (max(v,-Inf) > network.size(x))) 
    stop("Illegal vertex in remove.activity argument v.\n")
    
  # if deleting all vertices, can use network's delete.vertex.attribute
  # function, otherwise need to manually remove activity matrices.
  if (length(v) == network.size(x)) {
    delete.vertex.attribute(x, "active")
  } else {
    leave.active = setdiff(seq_along(x$val), v)
    leave.active = leave.active[!sapply(x$val[leave.active], is.null)]  #Filter out non-vertices
    left.activity = lapply(x$val[leave.active], "[[", "active")
    delete.vertex.attribute(x, "active")
    set.vertex.attribute(x, "active", left.activity, leave.active)
  }

  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)  
}

  
#--------------------------------------------------------------
# this is a helper function to insert a single valid spell.
# valid means that (i) onset <= terminus and (ii) onset != Inf
# and (iii) terminus != -Inf
#
# @param
#    spells  : the 2x(number of spells) matrix of current
#              spells
#    onset   : the onset time of the spell to be inserted;
#              default=-Inf
#    terminus: the terminus time of the spell to be inserted;
#              default=Inf
#
# @return:
#    the updated spells
#------------------------------------------------------------------
insert.spell<-function(spells, onset=-Inf, terminus=Inf){
  # forget all the below, do it in C
  return(.Call(InsertSpell_R, spells, onset, terminus, FALSE));

  if (is.null(spells) || spells[1,1]== Inf || spells[1,2]==-Inf)
    new.spells <- matrix(c(onset, terminus), 1,2)
  else {
    # case where no work is needed
    afton<-which(onset>=spells[,1])
    befon<-which(terminus<=spells[,2])
    if(length(afton)*length(befon)!=0 && max(afton)==min(befon)){
      afton = max(afton)
      if(!(onset==terminus && terminus==spells[afton,2] &&
           spells[afton,1]<spells[afton,2]))
        return(spells)
    }
    # time for work
    ns = NROW(spells)
    if(onset>spells[ns,1]) {  
      spell.row = ns+1  # row that spell will go in
    } else {
      spell.row = min(which(onset<=spells[,1]))
    }
    # spell row adjustments (continuations)
    if (spell.row > 1 && onset<=spells[spell.row-1,2] &&
        ifelse(terminus!=onset,T,terminus!=spells[spell.row-1,2])) {
      spell.row = spell.row-1
      onset = spells[spell.row,1]
    }
    if(terminus>=spells[ns, 2]) 
      retain.row = ns+1   # rows that are retained
    else
      retain.row = min(which(terminus<spells[,2]))  
    # retain row adjustments (continuations and instaneous pts)
    if(retain.row <= ns && terminus>=spells[retain.row,1]) {
      terminus=spells[retain.row,2] 
      retain.row=retain.row+1
    }else if(retain.row >1 && terminus==spells[retain.row-1,1] &&
             spells[retain.row-1,1]==spells[retain.row-1,2]){
      retain.row=retain.row-1
    }
    new.spells = matrix(c(onset, terminus),1,2)
    if(spell.row > 1)
      new.spells = rbind(spells[1:(spell.row-1),], new.spells)
    if(retain.row <= ns)
      new.spells = rbind(new.spells, spells[retain.row:ns,])
  }
  new.spells
}

#--------------------------------------------------------------
# this is a helper function to delete a single valid spell.
# valid means that (i) onset <= terminus and (ii) onset != Inf
# and (iii) terminus != -Inf
#
# @param
#    spells  : the 2x(number of spells) matrix of current
#              spells, assumed not null
#    onset   : the onset time of the spell to be deleted;
#              default=-Inf
#    terminus: the terminus time of the spell to be deleted;
#              default=Inf
#
# @return:
#    the updated spells
#------------------------------------------------------------------
delete.spell<-function(spells, onset=-Inf, terminus=Inf){
  # case where no work is needed
  ns = NROW(spells)
  if(onset > spells[ns,2] || (onset==spells[ns,2] && onset!=spells[ns,1]))
    return(spells)
  if(terminus < spells[1,1] || (terminus==spells[1,1] && onset!=terminus))
    return(spells)
  afton<-which(onset>=spells[,2])
  befon<-which(terminus<=spells[,1])
  if(length(afton)*length(befon)!=0 && max(afton)==min(befon)){
    afton = max(afton)
    if(!(spells[afton,1]==spells[afton,2] && onset==spells[afton,1]))
      return(spells)
  }
  # deactivation of points (disallowed for points in intervals)
  if (onset==terminus){
    pton = which(onset==spells[,1])
    if(length(pton)>0 && spells[max(pton),1]==spells[max(pton),2]){
      if(ns==1)
        return(matrix(c(Inf,Inf),1,2))
      else
        return(spells[-max(pton),,drop=FALSE])
    }else{
      return(spells)
    }
  }
  # deactivation of intervals
  if(onset<=spells[1,1])
    erow = 0   # the row number of the earlier rows to save
  else
    erow = max(which(onset>spells[,1]))
  if(terminus>spells[ns,2] || (terminus==spells[ns,2] &&
                               spells[ns,2]!=spells[ns,1])){
    lrow = ns+1   # the row number of the later rows to save
  } else if (terminus==spells[ns,2] & 
               spells[ns,2]==spells[ns,1]) {
    # when last spell is a matching point interval, keep it
    lrow <- ns
  } else {
    lrow = min(which(terminus<spells[,2]))
    if(lrow>1 && spells[lrow-1,2]==terminus && spells[lrow-1,1]==terminus)
      lrow = lrow-1
  }
  # erow and lrow adjustments (truncations)
  splitspell=NULL
  if(lrow==erow)
    splitspell = matrix(spells[lrow,],1,2)
  if(erow!=0 && onset < spells[erow,2])
    spells[erow,2] = onset
  if(lrow!=(ns+1) && terminus > spells[lrow,1]){
    if(lrow==erow){  # divided activation interval
      splitspell[1,1]=terminus
      lrow=lrow+1
    }else{
      spells[lrow,1] = terminus
    }
  }
  if(erow==0 && lrow==(ns+1))
    matrix(c(Inf,Inf),1,2)
  else if(erow==0)
    spells[(lrow:ns),,drop=FALSE]
  else if(lrow==(ns+1) && !is.null(splitspell))
    rbind(spells[(1:erow),,drop=FALSE], splitspell)
  else if(lrow==(ns+1))
    spells[(1:erow),,drop=FALSE]
  else if (erow+1==lrow && is.null(splitspell))
    spells
  else if (erow+1==lrow)
    rbind(spells[(1:erow),,drop=FALSE],
          splitspell,
          spells[(lrow:ns),,drop=FALSE])
  else
    spells[-((erow+1):(lrow-1)),,drop=FALSE]
}

# function to delete 'null' (Inf,Inf) spells from a lsit of spells
# and replace them with null
.removeNullSpells <- function(x){
  return(lapply(x,function(x){
    if(!is.null(x) && !is.na(x) && x[1,1]==Inf && x[1,2]==Inf){ return(NULL)} else { return(x) }
  }))
         
}
