#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2013 the statnet development team
######################################################################

#some basic functions to be used as a starting point for extending networkDynamic copied over from ndtv
# these may not fully meet spec
############ SOME UTILITY METHODS FOR WORKING WITH SPELLS, ATTRIBUTES, ETC #####



activate.network.attribute <- function (x, prefix, value, onset=NULL, terminus=NULL, length=NULL, at=NULL, dynamic.only=FALSE){
  
  #check that argument is a network
  if(!is.network(x)){
    stop("activate.network.attribute requires that the first argument be a network")
  }
  
  if(!is.character(prefix)){
    stop("prefixes must be character strings in activate.network.attribute.\n")
  }
  if(length(prefix) > 1){
    warning("Only the first element of prefix will be used.\n")
  }
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Activation times must be a numeric vector in activate.network.attribute.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in activate.network.attribute.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in activate.network.attribute.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in activate.network.attribute.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.logical(dynamic.only)){
    stop("dynamic.only flag must be a logical in activate.vertex.attribute.\n")
  }
  xn <- substitute(x) #this stuff is for modifying network inplace
  
  # figure out onset and terminus from at and length if necessary
  if(!is.null(at)) {
    onset <- terminus <- at
  } else if (!is.null(onset)) {

    if (!is.null(length))
      terminus <- onset + length
  } else {
    if (is.null(terminus)) {
      onset <- -Inf
      terminus <- Inf
    } else {
      onset <- terminus - length
    }
  }
  if(onset>terminus){
    stop("Onset times must precede terminus times in activate.network.attributes\n")
  }
  attrname = paste(prefix,"active",sep=".")
  
  # optionally replace non-active version of attribute with the same name
  if(!attrname%in%list.network.attributes(x)){
    if(prefix%in%list.network.attributes(x) & !dynamic.only){
      #delete the old attribute, it will effectively be replaced with the new one
      delete.network.attribute(x,prefix)
    }
  }
  
  
  timed <- get.network.attribute(x, attrname,unlist=FALSE); 
  if (!is.null(timed)){
    timed <- insertSpellAndVal(timed[[1]],timed[[2]],value,onset,terminus)
  } else {                                  #add the new dynamic attribute
    timed <- list(list(value),matrix(c(onset,terminus),nrow=1,ncol=2))
  }
  x <- set.network.attribute(x,attrname,value=timed)
  set.nD.class(x)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)	
}

#Only returns first value, since aggregation rules are not specified
get.network.attribute.active <- function(x, prefix, onset=NULL, terminus=NULL,length=NULL, at=NULL, rule = c("any", "all","earliest","latest"), dynamic.only=FALSE, return.tea=FALSE,unlist=TRUE){
  # checks for proper inputs
  if(!is.network(x)) 
    stop("get.network.attribute.active requires an argument of class network.\n")
  if(!is.null(at)) {
    if( !is.numeric(at))
      stop("'at' argument must be  numeric in get.network.attribute.active.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    if(length(at)!=1){
      stop("'at' time argument must have length one in get.network.attribute.active")
    }
  } else {
    if(!is.null(onset) &&  !is.numeric(onset))
      stop("Onset time must be  numeric  in get.network.attribute.active\n")
    if(!is.null(onset) &&  length(onset)!=1){
      stop("Onset time argument must have length one get.network.attribute.active\n")
    }
    if(!is.null(terminus) && !is.numeric(terminus))
      stop("Terminus times must be a numeric vector get.network.attribute.active.\n")
    if(!is.null(terminus) &&  length(terminus)!=1){
      stop("Terminus time argument must have length one get.network.attribute.active\n")
    }
    if(!is.null(length) && (!is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric in get.network.attribute.active\n")
    if(!is.null(length) &&  length(length)!=1){
      stop("Length time argument must have length one get.network.attribute.active\n")
    }
    
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(is.null(terminus) || is.null(length))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  
  if(!is.null(at)) {
    onset <- terminus <- at
  } else if (!is.null(onset)) {
    
    if (!is.null(length))
      terminus <- onset + length
  } else {
    if (is.null(terminus)) {
      onset <- -Inf
      terminus <- Inf
    } else {
      onset <- terminus - length
    }
  }
  
  if(onset>terminus)
    stop("Onset times must precede terminus times\n")
  
  rule<-match.arg(rule)
  
  if(return.tea){
    unlist=FALSE
    dynamic.only=TRUE
  }
  
  attrname = paste(prefix,"active",sep=".")
  
  # optionally replace non-active version of attribute with the same name
  if(!attrname%in%list.network.attributes(x)){
    if(prefix%in%list.network.attributes(x) & !dynamic.only){
      val <-get.network.attribute(x,prefix,unlist=unlist)
    } else {
      val <-NULL
    }
  } else {
    attribute <- get.network.attribute(x,attrname,unlist=FALSE)
    # find intersecting spells and corresponding values
    val <-getValsForSpells(list(attribute), onset=onset, terminus=terminus,
                          rule=rule, return.tea=return.tea)[[1]]
    # check for multiple values and return.tea state
    if (!return.tea){
      if(length(val)>1){
        warning("Multiple attribute values matched query spell for network attribute, only earliest value used\n")
        val <-val[[1]]
      }
      # TODO: this is probably where we would apply any other merge rules
      # now that there is only one value, need to remove a level of nesting.
    val <-val[[1]]
   }  
   if (unlist){
     val <- unlist(val)
   }
  }
  return(val)
}


get.edge.attribute.active<-function(x, prefix, onset=NULL, terminus=NULL,length=NULL, at=NULL, rule = c("any", "all","earliest","latest"), active.default = TRUE, dynamic.only=FALSE, require.active=FALSE,return.tea=FALSE,unlist=TRUE){
  if (!is.network(x)){
    stop("The first argument for get.edge.attribute.active must be a network object. (This version does not yet support using a list of edges)")
  }
  get.edge.value.active(x,prefix,onset,terminus,length,at,rule,active.default,dynamic.only,require.active,return.tea,unlist)
}

get.edge.value.active <- function(x, prefix, onset=NULL, terminus=NULL,length=NULL, at=NULL, rule = c("any", "all","earliest","latest"), active.default = TRUE, dynamic.only=FALSE, require.active=FALSE,return.tea=FALSE,unlist=TRUE){
  # validate inputs
  
  #check for common error of calling get.edge.attribute instead of get.edge.value
  if (!is.network(x)){
    stop("The first argument for get.edge.value.active must be a network object")
  }
  if (missing(prefix)){
    stop("The 'prefix' argument to get.value.attribute.active is missing, unable to determine the name of the attribute to fetch")
  }
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Singular time points given by the 'at' argument must be a numeric vector in get.edge.attribute.active\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in get.edge.value.active\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in get.edge.value.active\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in get.edge.value.active\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(is.null(terminus) || is.null(length))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  
  # process inputs
  rule<-match.arg(rule)
  if(!is.null(at)) {
    onset <- terminus <- rep(at, length=network.edgecount(x))
  } else if (!is.null(onset)) {
    onset <- rep(onset, length=network.edgecount(x))
    if(!is.null(terminus))
      terminus <- rep(terminus, length=network.edgecount(x))
    else if (!is.null(length))
      terminus <- onset + rep(length, length=network.edgecount(x))
  } else {
    terminus <- rep(terminus, length=length(el))
    onset <- terminus - rep(length, length=network.edgecount(x))
  }
  if(any(onset>terminus))
    stop("Onset times must precede terminus times\n")
  
  if(return.tea){
    unlist=FALSE
    dynamic.only=TRUE
  }
  
  attrname <- paste(prefix,"active",sep=".")
  #implement fall-through to non-tea
  if (!attrname%in%list.edge.attributes(x) & !dynamic.only){
    #check for non active version
    vals <- get.edge.value(x,prefix,unlist=FALSE)
  } else {
    vals<- as.list(rep(NA,network.edgecount(x)))
    attributes <- get.edge.value(x,attrname,unlist=FALSE)
  # locate potentially missing activity attributes
    valid.edges <- !sapply(attributes,is.null)
    # find intersecting spells and corresponding values
    vals[valid.edges] <-getValsForSpells(attributes[valid.edges], onset=onset, terminus=terminus,
                            rule=rule, return.tea=return.tea)
    # check for multiple values and return.tea state
    if (!return.tea){
      if(any(sapply(vals,length)>1)){
        warning(paste("Multiple attribute values matched query spell for attribute",attrname,"on some edges. Only earliest value used\n"))
        extravals <- which(sapply(vals[valid.edges],length)>1)
        for (index in extravals){
          vals[[index]]<-vals[[index]][[1]]
        }
      }
      # now that there is only one value, need to remove a level of nesting.
      vals[valid.edges] <-lapply(vals[valid.edges],"[[",1)
    }
  }
  # remove attributes for non-active edges depending on require.active
  if (require.active){
    vals[which(!is.active(x,onset=onset,terminus=terminus,e=seq_len(network.edgecount(x)),rule=rule,active.default=active.default))] <- NA
  }
  
  if (unlist){vals <-unlist(vals)}
  return(vals)
}

activate.vertex.attribute <- function (x, prefix, value, onset=NULL, terminus=NULL,length=NULL,at=NULL, v=seq_len(network.size(x)), dynamic.only=FALSE){
  #possibly duplicate list for when setting multiple verticies
  #check that argument is a network
  if(!is.network(x)){
    stop("activate.vertex.attribute requires that the first argument be a network")
  }
  if(!is.character(prefix)){
    stop("prefixes must be character strings in activate.vertex.attribute.\n")
  }
  if(length(prefix) > 1){
    warning("Only the first element of prefix will be used.\n")
  }
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Activation times must be a numeric vector in activate.vertex.attribute.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in activate.vertex.attribute.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in activate.vertex.attribute.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in activate.vertex.attribute.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.vector(v) || !is.numeric(v)){
    stop("Vertex ID's, v, must be a numeric vector in activate.vertex.attribute.\n")
  }
  if(!is.logical(dynamic.only)){
    stop("dynamic.only flag must be a logical in activate.vertex.attribute.\n")
  }
  if(length(v)==0){
  # no vertices specified so don't do anything  
   return(invisible(x))
  }
  if((min(v,Inf) < 1) || (max(v,-Inf) > network.size(x))){  # combinatorial breakage:
    stop("Illegal vertex id specified in v in activate.edge.attribute.\n")
  }
  xn <- substitute(x) #this stuff is for modifying network inplace
  
  # convert at and lenth options to onsets, and possibly replicate
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
  if(any(onset>terminus)){
    stop("Onset times must precede terminus times in activate.vertex.attribute.\n")
  }
  
  # TODO: also need to handle list case?
  if (length(value) != length(v)){
    if (!is.null(value)){
      value <- rep(value, length = length(v))
    }
  }
  attrname = paste(prefix,"active",sep=".")
  
  # optionally replace non-active version of attribute with the same name
  if(!attrname%in%list.vertex.attributes(x)){
    if(prefix%in%list.vertex.attributes(x) & !dynamic.only){
      #delete the old attribute, it will effectively be replaced with the new one
      delete.vertex.attribute(x,prefix)
    }
  }
  
  timedlist <- get.vertex.attribute(x, attrname,unlist=FALSE);
  # possibly get rid of values we not using
  timedlist <- timedlist[v]
  
  #have to loop instead of just checking if attribute exists because it can exist for some nodes and not others

  timedlist <-lapply(seq_len(length(v)),function(n){
    if (is.null(timedlist[[n]]) | (length(timedlist[[n]])==1 && is.na(timedlist[[n]]))){  #have to check if it is null or na, but na can be long
      return(list(list(value[[n]]),matrix(c(onset[n],terminus[n]),nrow=1,ncol=2)))
    } else {
      return(insertSpellAndVal(timedlist[[n]][[1]],timedlist[[n]][[2]],value[[n]],onset[[n]],terminus[[n]]))
      
    }
  })
  
  x <- set.vertex.attribute(x,attrname,value=timedlist,v=v)
  set.nD.class(x)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)	
}


activate.edge.value<-function(x, prefix, value, onset=NULL, terminus=NULL,length=NULL,at=NULL, 
                                   e=seq_along(x$mel), dynamic.only=FALSE){
  
  if (!is.network(x)) 
    stop("activate.edge.value requires an argument of class network.\n")
  if (is.hyper(x)) 
    stop("Hypergraphs not currently supported in activate.edge.value.\n")
  n <- network.size(x)
  if (!is.matrix(value)) {
    if (is.vector(value)) 
      value <- matrix(rep(value, length = n * n), n, n)
    else value <- matrix(value, n, n)
  }
  if ((min(e) < 1) | (max(e) > length(x$mel))) 
    stop("Illegal edge in activate.edge.value.\n")
  xn <- substitute(x)
  
  # remove e corresponding to NULL edge values (deleted edges)
  e<-e[!sapply(x$mel[e],is.null)]
 
  # get the list of values corresponding to edges that exist, limited by e
  evals<- lapply(x$mel[e],function(edge){
    return(value[edge$outl,edge$inl])
  })
  
  x<-activate.edge.attribute(x=x,prefix=prefix,value=evals, onset=onset,terminus=terminus,length=length,at=at,e=e,dynamic.only=dynamic.only)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)
  
}

#this may be too slow if edge allready has spells for the attribute
# also stores its vals in a vector, needs to be a list
activate.edge.attribute <-function(x, prefix, value, onset=NULL, terminus=NULL,length=NULL,at=NULL, 
	e=seq_along(x$mel), dynamic.only=FALSE){
  #check that argument is a network
  if(!is.network(x)){
    stop("activate.edge.attribute requires that the first argument be a network")
  }
  #check that it has edges
  if(network.edgecount(x)==0){
    warning("edge attributes cannot be activated because network does not contain any edges")
  }
  
  if(!is.character(prefix)){
    stop("prefixes must be character strings in activate.edge.attribute.\n")
  }
  if(length(prefix) > 1){
    warning("Only the first element of prefix will be used.\n")
  }
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Activation times must be a numeric vector in activate.edge.attribute.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in activate.edge.attribute.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in activate.edge.attribute.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in activate.edge.attribute.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.vector(e) || !is.numeric(e)){
    stop("Edge ID's, e, must be a numeric vector in activate.edge.attribute.\n")
  }
  if(!is.logical(dynamic.only)){
    stop("dynamic.only flag must be a logical in activate.edge.attribute.\n")
  }
  if (length(e)>0){
    if((min(e,Inf) < 1) || (max(e,-Inf) > length(x$mel))){  # combinatorial breakage:
      stop("Illegal edge id specified in e in activate.edge.attribute.\n")
    }
  }
  
  xn <- substitute(x) #this stuff is for modifying network inplace
  
  # possibly convert at and length or replicated onsets
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
  if(any(onset>terminus)){
    stop("Onset times must precede terminus times in activate.edge.attribute\n")
  }
  # replicate values if they were not enough per edge
  if (length(value) != length(e)){
    value <- rep(value, length = length(e))
  }
  
  #Filter outvalues corresponding to non-edges caused by edge deletion
  value <-value[!sapply(x$mel[e], is.null)]
  e <- e[!sapply(x$mel[e], is.null)]  
  
	attrname <- paste(prefix,"active",sep=".")	
  
  # optionally replace non-active version of attribute with the same name
  knownAttrs<-list.edge.attributes(x)
  if(!attrname%in%knownAttrs){
    if(prefix%in%knownAttrs & !dynamic.only){
      #delete the old attribute, it will effectively be replaced with the new one
      delete.edge.attribute(x,prefix)
    }
  }
  
  
  timedlist <- get.edge.attribute(x$mel, attrname,unlist=FALSE);
  timedlist <- lapply(seq_len(length(e)),function(n){
    timed <- timedlist[[e[n]]];
    if(is.null(timed) || (length(timed)==1 && is.na(timed))){
      # create a new TEA attribute
      return(list(list(value[[n]]),matrix(c(onset[n],terminus[n]),nrow=1,ncol=2)))
    } else {
      # append to existing TEA attribute
      return(insertSpellAndVal(timed[[1]],timed[[2]],value[[n]],onset[[n]],terminus[[n]]))
    }
  })
  set.edge.attribute(x,attrname,timedlist,e=e)
	set.nD.class(x)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
	invisible(x)
}


get.vertex.attribute.active <- function(x, prefix, onset=NULL, terminus=NULL,length=NULL, at=NULL, rule = c("any", "all","earliest","latest"), na.omit = FALSE, null.na = TRUE, active.default = TRUE, dynamic.only=FALSE, require.active=FALSE,return.tea=FALSE,unlist=TRUE){
  
  # validate inputs
  if(!is.network(x)) 
    stop("get.vertex.attribute.active requires an argument of class network.\n")
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Singular time points given by the 'at' argument must be a numeric vector in get.vertex.attribute.active\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in get.vertex.attribute.active\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in get.vertex.attribute.active\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in get.vertex.attribute.active\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(is.null(terminus) || is.null(length))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  
  # process inputs
  rule<-match.arg(rule)
  if(!is.null(at)) {
    onset <- terminus <- rep(at, length=network.size(x))
  } else if (!is.null(onset)) {
    onset <- rep(onset, length=network.size(x))
    if(!is.null(terminus))
      terminus <- rep(terminus, length=network.size(x))
    else if (!is.null(length))
      terminus <- onset + rep(length, length=network.size(x))
  } else {
    terminus <- rep(terminus, length=network.size(x))
    onset <- terminus - rep(length, length=network.size(x))
  }
  if(any(onset>terminus))
    stop("Onset times must precede terminus times\n")
  
  if(return.tea){
    unlist=FALSE
    dynamic.only=TRUE
  }
  
  attrname <- paste(prefix,"active",sep=".")
  if (!attrname%in%list.vertex.attributes(x) & !dynamic.only){
    #check for non active version
    vals <- get.vertex.attribute(x,prefix,na.omit=na.omit,null.na=TRUE,unlist=FALSE)
  } else {
    attributes <- get.vertex.attribute(x,attrname,na.omit=na.omit,null.na=TRUE,unlist=FALSE)
    # find intersecting spells and corresponding values
    vals <-getValsForSpells(attributes, onset=onset, terminus=terminus, 
                                      rule=rule, return.tea=return.tea)
    # check for multiple values and return.tea state
    if (!return.tea){
      extravals <- which(sapply(vals,length)>1)
      if(length(extravals)>0){
        warning(paste("Multiple attribute values matched query spell for attribute",attrname,"on some vertices. Only earliest value used\n"))
        for (index in extravals){
          vals[[index]]<-vals[[index]][[1]]
        }
      }
      # TODO: this is probably where we would apply any other merge rules
      # now that there is only one value, need to remove a level of nesting.
      vals <-lapply(vals,"[[",1)
    }
  }
  # remove attributes for non-active vertices depending on require.active
  if (require.active){
    vals[which(!is.active(x,onset=onset,terminus=terminus,v=seq_len(network.size(x)),rule=rule,active.default=active.default))] <- NA
  }
  #  handle null.na
  # this is awfully awkward, but not sure anyone uses anyways
  if (!null.na){
    possiblyNull <- get.vertex.attribute(x,attrname,na.omit=na.omit,null.na=FALSE,unlist=FALSE)
    vals[sapply(possiblyNull, is.null)] <- list(NULL)
  }
  # todo: handle na.omit and null.na here!
  if (unlist){vals <-unlist(vals)}
  return(vals)
}


#function to check if two spells intersect
spells.overlap <-function(s1,s2){
  if (length(s1)!=2 | length(s2)!=2){
		stop("Each spell must be a vector of length 2")
	}
  # check for null spell
  if ((s1[1]==Inf & s1[2]==Inf) | (s2[1]==Inf & s2[2]==Inf)){
    return(FALSE)
  }
  # check for equality (point -point and spell spell)
  if (all(s1==s2)){
    return(TRUE) 
  }
  # point-spell comparisons
  if (s1[1]==s1[2]){
    if (s1[1] >= s2[1] & s1[1] < s2[2]){
      return(TRUE)
    }
  } else if (s2[1]==s2[2]) {
    if (s2[1] >= s1[1] & s2[1] < s1[2]){
      return(TRUE)
    }
  }
  
  # interval comparisons
  #start one < end two and end one > start 2
  if (s1[1]<s2[2] & s1[2] > s2[1]){
    return(TRUE)
  }
  
	return(FALSE)  # non-overlapping
	
}

#search an array of spells to see if any intersect
spells.hit<-function(needle,haystack){
  if(length(needle)!=2){
    stop("search spell (needle) must have exactly two elements")
  }
  if(ncol(haystack)!=2){
    stop("target  spell matrix (haystack) must have exactly two columns")
  }
  if(nrow(haystack)<1){
    stop("target  spell matrix (haystack) must have at least one row")
  }
  for (s in 1:nrow(haystack)){
    if(spells.overlap(needle,haystack[s,])){
      return(s)
    }
  }
  return(-1)
}



search.spell<-function(needle,haystack){
  if(length(needle)!=2){
    stop("search spell (needle) must have exactly two elements")
  }
  if(needle[1]>needle[2]){
    stop("search spell (needle) must have onset <= terminus")
  }
  if(ncol(haystack)!=2){
    stop("target  spell matrix (haystack) must have exactly two columns")
  }
  if(nrow(haystack)<1){
    stop("target  spell matrix (haystack) must have at least one row")
  }
  indices <- which(sapply(1:nrow(haystack),function(s){ spells.overlap(needle,haystack[s,])}))
  return(indices)
}



# Ayn wrote this, which is why it is so nicely documented
#--------------------------------------------------------------
# this is a helper function to insert a single spell
#
# @param
#    cur.vals  : a list of current values
#    cur.spells: the 2x(number of spells) matrix of current
#                spells
#    val       : the value of the spell to be inserted
#    onset     : the onset time of the spell to be inserted;
#                default=-Inf
#    terminus  : the terminus time of the spell to be inserted;
#                default=Inf
#
# @return:
#    the list of two elements, 1) the updated values, 2) the
#    updated spells
#------------------------------------------------------------------

#for this version, try no testing, just do operation
fastInsertSpellAndVal<-function(cur.vals, cur.spells, val, onset=-Inf,
         terminus=Inf){
  list(
    c(cur.vals, list(val)),
    matrix(c(cur.spells,onset,terminus),ncol=2,byrow=TRUE)
  )
}


temp.insertSpellAndVal<-function(cur.vals, cur.spells, val, onset=-Inf,
                       terminus=Inf){
     if (is.null(cur.spells) ||
         (is.infinite(cur.spells[1,1]) & cur.spells[1,1] > 0)) {
       new.values <- list(val)
       new.spells <- matrix(c(onset, terminus), 1,2)
     } else {
        double.entry=FALSE
        ns <- NROW(cur.spells)
        # find the row that this spell will go into
      	if(onset>cur.spells[ns,1]) {
      	  spell.row <- ns+1
      	} else {
      	  spell.row <- min(which(onset<=cur.spells[,1]))
      	}
      	# if the onset interrupts/continues an existing spell... 
      	if (spell.row > 1 && onset<=cur.spells[spell.row-1,2]) {
                if(identical(val, cur.vals[[spell.row-1]])) {
      	     spell.row <- spell.row-1
      	     onset <- cur.spells[spell.row,1]   # back up onset to that of interrupted spell
      	   } else {
                   if(terminus < cur.spells[spell.row-1,2]){  # this spells splits the interrupted spell in 2
                     double.entry<-TRUE
                     val2 <- cur.vals[[spell.row-1]]
                     spell2 <-matrix(c(terminus, cur.spells[spell.row-1,2]),1,2)
                   }
                   cur.spells[spell.row-1,2] <- onset   # truncate interrupted spell
      	   }
      	}
      	# find the minimum spell that is retained (vs. spells that are overlapped/overwritten)
      	if(terminus>=cur.spells[ns, 2]){
      	  retain.row <- ns+1
      	} else {
      	  retain.row <- min(which(terminus<cur.spells[,2]))
      	}
      	# if the terminus interrupts/continues an existing spell...
      	if(retain.row <= ns && terminus>=cur.spells[retain.row,1]) {
      	  if(identical(val,cur.vals[[retain.row]])) {
      	    terminus <- cur.spells[retain.row,2]  # forward terminus to that of interrupted spell
      	    retain.row <- retain.row+1
      	  } else {
      	    cur.spells[retain.row,1] <- terminus  # partial overwrite of interrupted spell
      	  }
      	}
      	# construct new spell matrix and value vector
      	#new.spells <- matrix(c(onset, terminus),1,2)
      	#new.values <- list(val)
              if(double.entry){
                new.spells <- rbind(new.spells, spell2)
                new.values[[2]] <- val2
              }
      	if(spell.row > 1){
      	  #new.spells <- rbind(cur.spells[1:(spell.row-1),], new.spells)
          new.spells<-matrix(c(as.numeric(cur.spells[1:(spell.row-1),]),onset,terminus),nrow=spell.row,ncol=2,byrow=TRUE)
      	  new.values <- c(cur.vals[1:(spell.row-1)], list(val))
      	}
      	if(retain.row <= ns) {
      	  new.spells <- rbind(new.spells, cur.spells[retain.row:ns,])
      	  new.values <- c(new.values, cur.vals[retain.row:ns])
      	}
    }
    list(new.values, new.spells)
}

insertSpellAndVal<-function(cur.vals, cur.spells, val, onset=-Inf,
                            terminus=Inf){
  if (is.null(cur.spells) ||
        (is.infinite(cur.spells[1,1]) & cur.spells[1,1] > 0)) {
    new.values <- list(val)
    new.spells <- matrix(c(onset, terminus), 1,2)
  } else {
    double.entry=FALSE
    ns <- NROW(cur.spells)
    # find the row that this spell will go into
    if(onset>cur.spells[ns,1]) {
      spell.row <- ns+1
    } else {
      spell.row <- min(which(onset<=cur.spells[,1]))
    }
    # if the onset interrupts/continues an existing spell... 
    if (spell.row > 1 && onset<=cur.spells[spell.row-1,2]) {
      if(identical(val, cur.vals[[spell.row-1]])) {
        spell.row <- spell.row-1
        onset <- cur.spells[spell.row,1]   # back up onset to that of interrupted spell
      } else {
        if(terminus < cur.spells[spell.row-1,2]){  # this spells splits the interrupted spell in 2
          double.entry<-TRUE
          val2 <- cur.vals[[spell.row-1]]
          spell2 <-matrix(c(terminus, cur.spells[spell.row-1,2]),1,2)
        }
        cur.spells[spell.row-1,2] <- onset   # truncate interrupted spell
      }
    }
    # find the minimum spell that is retained (vs. spells that are overlapped/overwritten)
    if(terminus>=cur.spells[ns, 2]){
      retain.row <- ns+1
    } else {
      retain.row <- min(which(terminus<cur.spells[,2]))
    }
    # if the terminus interrupts/continues an existing spell...
    if(retain.row <= ns && terminus>=cur.spells[retain.row,1]) {
      if(identical(val,cur.vals[[retain.row]])) {
        terminus <- cur.spells[retain.row,2]  # forward terminus to that of interrupted spell
        retain.row <- retain.row+1
      } else {
        cur.spells[retain.row,1] <- terminus  # partial overwrite of interrupted spell
      }
    }
    # construct new spell matrix and value vector
    new.spells <- matrix(c(onset, terminus),1,2)
    new.values <- list(val)
    if(double.entry){
      new.spells <- rbind(new.spells, spell2)
      new.values[[2]] <- val2
    }
    if(spell.row > 1){
      new.spells <- rbind(cur.spells[1:(spell.row-1),], new.spells)
      new.values <- c(cur.vals[1:(spell.row-1)], new.values)
    }
    if(retain.row <= ns) {
      new.spells <- rbind(new.spells, cur.spells[retain.row:ns,])
      new.values <- c(new.values, cur.vals[retain.row:ns])
    }
  }
  list(new.values, new.spells)
}

# define functions to use as match rules. These are interval comparison functions 
# x = a spell [onset,terminus]
# y = onset to be compared for match
# z = terminus to be compared for match
match.rule.all <-function(x,y,z){y >= x[1] && z <= x[2]}
match.rule.any <-function(x,y,z){(y >= x[1] && y <  x[2]) ||
                              (z >  x[1] && z <= x[2]) ||
                (y <= x[1] && z >= x[2] && !(x[2]==x[1] && z==x[2]))}


# internal function to get a set of values corresponding to a query spell, uing the TEA value and spell list components
# active is a list of TEA attributes corresponding to a vertex (or edge), each of which is a list where the first element is list of atribute values, and the second element 2nd is spell matrix
# onset, terminus, length and at if specified must have the same length as active
# returns a list with elements corresponding to active, which are NA or contain list of matching values
# TODO: maybe this could be re-written to be faster with the "findInterval" function?
# or just move the entire function into C
# or parallelize with package "parallel" and parLapply
getValsForSpells <- function(active,onset=NULL,terminus=NULL,length=NULL, at=NULL,
                             rule=c("any","all","earliest","latest"),return.tea=FALSE){
  
  # since this is not a user function assumes inputs were allready verified by the calling function

  # switch out the function we will use to compare spell intervals 
  # to either the 'any' or 'all' match rule
  if(match.arg(rule)=="all"){
    int.query.true <- match.rule.all
  } else {
    int.query.true <- match.rule.any
  }
  
  
  # create a vector default values of NA for any elements that may not be otherwise specified
  vals<-as.list(rep(NA,times=length(active)))
  
  # Skye: I tried this loop as lapply, didn't get any speedup, but that could be because it wasn't well written. 
  
  # for each network element (either vertex or edge), check the spell matrix for intersections with the query spell and return the corresponding values. 
  for(i in seq_along(active)){
    if(is.null(active[[i]])){ # is the attribute malformed?
      stop(paste("activity attribute missing for element",i))
      # TODO: should this be controled by the null.na value?
    } else if (length(active[[i]])==1 && is.na(active[[i]])){ # if no value was ever specified ...
      # don't do anything, we will return NA
    } else if (active[[i]][[2]][1,1]==-Inf && active[[i]][[2]][1,2]==Inf){ # if it is always active ...
      if (return.tea){ # should we return a TEA structure (for user to process) instead of just the value?
        vals[[i]]<-list(active[[i]][[1]][[1]],active[[i]][[2]][1,])
      } else {
        vals[[i]]<-active[[i]][[1]][[1]]  # just return the value
      }
    } else if (all(is.infinite(active[[i]][[2]]))) { # check allways inactive spell (Inf, Inf) The valid (-Inf,Inf) already dealt with in previous condition.
      vals[[i]]<-NA  # it was marked as never active, so set to NA
    } else if (terminus[i] == onset[i]){  # if the query is a point interval we know it can only match one value
      # construct a vector of indices of spells having an onset equal to the query onset
      befon<-which(onset[i]==active[[i]][[2]][,1]) # before onset
      if(length(befon)>0){  # did we find it?
        if (return.tea){ # should we return a TEA structure (for user to process) instead of just the value?
          vals[[i]]<-list(active[[i]][[1]][befon],active[[i]][[2]][befon,])
        } else {
          vals[[i]]<-active[[i]][[1]][befon]
        }
      } else {
        # construct a vector of indices of spells where the query onset is after (>=) spell onset
        afton <-which(onset[i]>=active[[i]][[2]][,1])
        # construct a vector of indices of spells where the query onset is before (<) spell terminus
        befterm<-which((active[[i]][[2]][,2]==Inf)|(onset[i]<active[[i]][[2]][,2]))
        # intersect the vector of indices to find the matching set
        splindex <-sort(intersect(afton,befterm))
        if (length(splindex)>0){
          if (return.tea){ # should we return a TEA structure (for user to process) instead of just the value?
            vals[[i]]<-list(active[[i]][[1]][splindex],active[[i]][[2]][splindex,])
          } else {
            vals[[i]]<-active[[i]][[1]][splindex]
          }
        }
      }
    } else {  # interval query
      # TODO: this seems to be slow..
      # apply the spell query to each row of the activity matrix to determine the indices of intersecting values
      # old code:
      #splindex <- which(apply(active[[i]][[2]], 1, int.query.true, onset[i], terminus[i]))
      # new code:
      # amazingly, the for loop below takes less time than apply above because we can loop backwards
      splindex<-numeric(0)
      foundSome<-FALSE
      for (r in nrow(active[[i]][[2]]):1){
        if(int.query.true(active[[i]][[2]][r,],onset[i],terminus[i])){
          splindex<-c(r,splindex)
          foundSome<-TRUE
        } else if (foundSome){ # since rows are ordered, stop searching if we've already found some
          break
        }
      }
      if (length(splindex)>0){  # if we found 1 or more values ...
        # evaluate the earliest or latest rule to decide which to return
        if(rule=='earliest'){
          splindex<-splindex[1] # choose the earliest value found
        } else if (rule=='latest'){
          splindex<-splindex[length(splindex)] # choose the latest value found
        }
        if (return.tea){ # should we return the TEA structure instead of just the value?
          vals[[i]]<-list(active[[i]][[1]][splindex],active[[i]][[2]][splindex,])
        } else {
          vals[[i]]<-active[[i]][[1]][splindex]
        }
      } 
    }
  }
  vals
}

#questions?
# 1. after deactivate all spells, should we leave the tea name with an empty list or delete the tea completely?
# 
##########################################################################
#deactivate.vertex.attribute
##########################################################################
deactivate.vertex.attribute <- function (x, prefix, onset=NULL, terminus=NULL,length=NULL,at=NULL, v=seq_len(network.size(x)), dynamic.only=FALSE){
  #possibly duplicate list for when setting multiple verticies
  #check that argument is a network
  if(!is.network(x)){
    stop("deactivate.vertex.attribute requires that the first argument be a network")
  }
  if(!is.character(prefix)){
    stop("prefixes must be character strings in deactivate.vertex.attribute. \n")
  }
  if(length(prefix) > 1){
    warning("Only the first element of prefix will be used.\n")
  }
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Activation times must be a numeric vector in deactivate.vertex.attribute. \n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in deactivate.vertex.attribute. \n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in deactivate.vertex.attribute. \n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in deactivate.vertex.attribute. \n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.vector(v) || !is.numeric(v)){
    stop("Vertex ID's, v, must be a numeric vector in deactivate.vertex.attribute. \n")
  }
  if(!is.logical(dynamic.only)){
    stop("dynamic.only flag must be a logical in deactivate.vertex.attribute.\n")
  }
  if((min(v,Inf) < 1) || (max(v,-Inf) > network.size(x))){  # combinatorial breakage:
    stop("Illegal vertex id specified in v in deactivate.vertex.attribute.\n")
  }
  xn <- substitute(x) #this stuff is for modifying network inplace
  
  # convert at and lenth options to onsets, and possibly replicate
  if(!is.null(at)) {
    onset <- terminus <- rep(at, length=network.size(x)) # modifed on 04/04 from legnth(v) to network.size(x))
  } else if (!is.null(onset)) {
    onset <- rep(onset, length=network.size(x))
    if(!is.null(terminus))
      terminus <- rep(terminus, length=network.size(x))
    else if (!is.null(length))
      terminus <- onset + rep(length, length=network.size(x))
  } else {
    if (is.null(terminus)) {
      onset <- rep(-Inf, length=network.size(x))
      terminus <- rep(Inf, length=network.size(x))
    } else {
      terminus <- rep(terminus, length=network.size(x))
      onset <- terminus - rep(length, length=network.size(x))
    }
  }
  if(any(onset>terminus)){
    stop("Onset times must precede terminus times in deactivate.vertex.attribute.\n")
  }
  attrname = paste(prefix,"active",sep=".")
  
  # optionally replace non-active version of attribute with the same name
  if(!attrname%in%list.vertex.attributes(x)){
    if(prefix%in%list.vertex.attributes(x) & !dynamic.only){
      #delete the old attribute, it will effectively be replaced with the new one
      delete.vertex.attribute(x,prefix)
    }
  }
  
  timedlist <- get.vertex.attribute(x, attrname,unlist=FALSE);
  # possibly get rid of values we not using
  
  # 4/1,fix problem for v is not entire network nodes.
  
#   if(is.null(timedlist)){
#     warning("intend to deactivate attribute on inactivate vertex")
#     return(character(0))}
#   v.act <- which(!unlist(lapply(timedlist,function(x)all(is.na(x)))))
#   if (!all(v%in%v.act)){
#     warning("intend to deactivate attribute on inactivate vertex")
#   v = intersect(v.act,v)}
#   
   timedlist <- timedlist[v]
  
  #have to loop instead of just checking if attribute exists because it can exist for some nodes and not others
  
  timedlist <-lapply(seq_len(length(timedlist)),function(n){
    if (!is.null(timedlist[[n]]) && !(length(timedlist[[n]])==1) && !is.na(timedlist[[n]])){ 
      return(
        deactive.spell.attribute(onset=onset[n],terminus=terminus[n],spell.mat=timedlist[[n]][[2]],val.list=timedlist[[n]][[1]])
      )}})
  
  x <- set.vertex.attribute(x,attrname,value=timedlist,v=v)
  set.nD.class(x)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)  
}


##########################################################################
#withIn: a function to select the changing spell
##########################################################################
withIn <- function(onset,terminus,spell.mat){
  apply(spell.mat,1,function(x) !(terminus<=x[1] | onset>=x[2]))
}

##########################################################################
#deactive.spell.attribute:  a internal function to output the spell matrix
#                           after deactivate
##########################################################################
deactive.spell.attribute <- function(onset,terminus,spell.mat,val.list){
  index.row <- withIn(onset,terminus,spell.mat)
  if(any(index.row)==F){ 
    return(list(val.list,spell.mat))
  } else{
    val.list.c <- lapply(which(index.row),function(x)val.list[[x]])
    add.row <- t(apply(spell.mat[index.row,,drop=FALSE],1,function(x){ 
      #     x<-spell.mat[index.row,,drop=FALSE]
      if(x[1]>=onset && x[2]>terminus) #deactivate at the beginning
        return(c(terminus,x[2]))
      else if(x[1]>=onset && x[2]<=terminus) # deactivate covers the spell
        return(c(NA,NA))
      else if(x[1]<onset && x[2]>terminus) # deactivate at the interval within the spell
        return(matrix(c(x[1],onset,terminus,x[2]),byrow=T,ncol=2))
      else if(x[1]<onset && x[2]<=terminus) # deactiate at the ending 
        return(c(x[1],onset))
    }))
    if(dim(add.row)[2]==4){
      # temp version to solve the matrix returns....may need to improve
      add.row <- matrix(add.row,ncol=2,byrow=F)
      val.list.c<-rep(val.list.c,2)}
    
    add.row <- na.omit(add.row)
    
    # for the spells being compeletely deactivated, NA,NA returned, and need to delete the corresponding val. 
    if(!is.null(attr(add.row,"na.action")))
      val.list.c <-val.list.c[!(1:length(val.list.c) %in% attr(add.row,"na.action"))] 
    
    spell.mat.new <- spell.mat[!index.row,]
    val.list.new <- val.list[!index.row]
    spell.mat.new <- rbind(spell.mat.new,add.row)
    val.list.new <-append(val.list.new,val.list.c)
    #   spell.mat.new <- spell.mat.new[order(val.new,spell.mat.new[,1]),]
    if (!(length(val.list.new) | length(spell.mat.new))){
      if(!length(val.list.new))
        val.list.new <- list("NA")
      if(!length(spell.mat.new))
        spell.mat.new <- matrix(c(Inf,Inf),ncol=2)} 
    else{
      val.list.new <- lapply(1:length(val.list.new),function(x)val.list.new[[order(spell.mat.new[,1],spell.mat.new[,2],decreasing=FALSE)[x]]])
      spell.mat.new <- unname(spell.mat.new[order(spell.mat.new[,1],spell.mat.new[,2],decreasing=FALSE),,drop=FALSE])}
    # to do: order the list, 
    
    return(list(val.list.new,spell.mat.new))
  }
}


##########################################################################
#deactivate.edge.attribute
##########################################################################

deactivate.edge.attribute <-function(x, prefix, onset=NULL, terminus=NULL,length=NULL,at=NULL,e=seq_along(x$mel), dynamic.only=FALSE){
  #check that argument is a network
  if(!is.network(x)){
    stop("deactivate.edge.attribute requires that the first argument be a network")
  }
  #check that it has edges
  if(network.edgecount(x)==0){
    warning("edge attributes cannot be deactivated because network does not contain any edges")
  }
  
  if(!is.character(prefix)){
    stop("prefixes must be character strings in deactivate.edge.attribute. \n")
  }
  if(length(prefix) > 1){
    warning("Only the first element of prefix will be used.\n")
  }
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Activation times must be a numeric vector in deactivate.edge.attribute.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in deactivate.edge.attribute.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in deactivate.edge.attribute.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in deactivate.edge.attribute.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.vector(e) || !is.numeric(e)){
    stop("Edge ID's, e, must be a numeric vector in deactivate.edge.attribute.\n")
  }
  if(!is.logical(dynamic.only)){
    stop("dynamic.only flag must be a logical in deactivate.edge.attribute.\n")
  }
  if (length(e)>0){
    if((min(e,Inf) < 1) || (max(e,-Inf) > length(x$mel))){  # combinatorial breakage:
      stop("Illegal edge id specified in e in deactivate.edge.attribute.\n")
    }
  }
  
  xn <- substitute(x) #this stuff is for modifying network inplace
  
  #Filter out non-edges caused by edge deletion
  e <- e[!sapply(x$mel[e], is.null)]  
  
  # possibly convert at and length or replicated onsets
  if(!is.null(at)) {
    onset <- terminus <- rep(at, length=length(seq_along(x$mel)))  # modifed on 04/04 from legnth(e) to length(seq_along(x$mel))
  } else if (!is.null(onset)) {
    onset <- rep(onset, length=length(seq_along(x$mel)))
    if(!is.null(terminus))
      terminus <- rep(terminus, length=length(seq_along(x$mel)))
    else if (!is.null(length))
      terminus <- onset + rep(length, length=length(seq_along(x$mel)))
  } else {
    if (is.null(terminus)) {
      onset <- rep(-Inf, length=length(seq_along(x$mel)))
      terminus <- rep(Inf, length=length(seq_along(x$mel)))
    } else {
      terminus <- rep(terminus, length=length(seq_along(x$mel)))
      onset <- terminus - rep(length, length=length(seq_along(x$mel)))
    }
  }
  if(any(onset>terminus)){
    stop("Onset times must precede terminus times in deactivate.edge.attribute.\n")
  }
  
  attrname = paste(prefix,"active",sep=".")  
  
  # optionally replace non-active version of attribute with the same name
  knownAttrs<-list.edge.attributes(x)
  if(!attrname%in%knownAttrs){
    if(prefix%in%knownAttrs & !dynamic.only){
      #delete the old attribute, it will effectively be replaced with the new one
      delete.edge.attribute(x,prefix)
    }
  }
  
  
  timedlist <- get.edge.attribute(x$mel,attrname,unlist=FALSE)
  
#   # 4/1,fix problem for e is not entire network edges.
#   
#   if(is.null(timedlist)){
#     warning("intend to deactivate attribute on inactivate edges")
#     return(character(0))}
#   e.act <- which(!unlist(lapply(timedlist,function(x)all(is.null(x)))))
#   if (!all(e%in%e.act)){
#     warning("intend to deactivate attribute on inactivate edge")
#     e = intersect(e.act,e)}
#   
  
#   timedlist <- timedlist[e]
  
  timedlist <-lapply(seq_len(length(e)),function(n){
    timed <- timedlist[[n]];
    if (!is.null(timed) && !(length(timed)==1) && !is.na(timed)){ 
      return(
        deactive.spell.attribute(onset=onset[n],terminus=terminus[n],spell.mat=timed[[2]],val.list=timed[[1]])
      )}})
  
  set.edge.attribute(x,attrname,timedlist,e=e)
  set.nD.class(x)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)
}


##########################################################################
# -----------deactivate.network.attribute----
##########################################################################


deactivate.network.attribute <- function (x, prefix, onset=NULL, terminus=NULL, length=NULL, at=NULL, dynamic.only=FALSE){
  
  if(!is.network(x)){
    stop("deactivate.network.attribute requires that the first argument be a network")
  }
  
  if(!is.character(prefix)){
    stop("prefixes must be character strings in deactivate.network.attribute.\n")
  }
  if(length(prefix) > 1){
    warning("Only the first element of prefix will be used.\n")
  }
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Activation times must be a numeric vector in deactivate.network.attribute.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in deactivate.network.attribute.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in deactivate.network.attribute.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in deactivate.network.attribute.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.logical(dynamic.only)){
    stop("dynamic.only flag must be a logical in deactivate.network.attribute.\n")
  }
  xn <- substitute(x) #this stuff is for modifying network inplace
  
  # figure out onset and terminus from at and length if necessary
  if(!is.null(at)) {
    onset <- terminus <- at
  } else if (!is.null(onset)) {
    
    if (!is.null(length))
      terminus <- onset + length
  } else {
    if (is.null(terminus)) {
      onset <- -Inf
      terminus <- Inf
    } else {
      onset <- terminus - length
    }
  }
  if(onset>terminus){
    stop("Onset times must precede terminus times in deactivate.network.attribute.\n")
  }
  attrname = paste(prefix,"active",sep=".")
  
  # optionally replace non-active version of attribute with the same name
  if(!attrname%in%list.network.attributes(x)){
    if(prefix%in%list.network.attributes(x) & !dynamic.only){
      #delete the old attribute, it will effectively be replaced with the new one
      delete.network.attribute(x,prefix)
    }
  }
  
  
  timed <- get.network.attribute(x, attrname,unlist=FALSE); 
  if (!is.null(timed) | !(length(timed)==1 && !is.na(timed))){ 
    timed <- 
      deactive.spell.attribute(onset=onset,terminus=terminus,spell.mat=timed[[2]],val.list=timed[[1]])
  }
  
  x <- set.network.attribute(x,attrname,value=timed)
  set.nD.class(x)
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)  
}



##########################################################################
#list.vertex.attributes.active
##########################################################################

list.vertex.attributes.active <-function(x, onset=NULL, terminus=NULL,length=NULL, at=NULL, na.omit=FALSE, rule = c("any", "all"), v=seq_len(network.size(x)), require.active=FALSE, active.default = TRUE, dynamic.only=FALSE){
  
  if(!is.network(x)){
    stop("list.vertex.attributes.active requires that the first argument be a network object")
  }
  
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Specified times must be a numeric vector in list.vertex.attributes.active\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in list.vertex.attributes.active.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in list.vertex.attributes.active.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in list.vertex.attributes.active.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  if(!is.vector(v) || !is.numeric(v)){
    stop("Vertex ID's, v, must be a numeric vector in list.vertex.attributes.active.\n")
  }
  if(!is.logical(dynamic.only)){
    stop("dynamic.only flag must be a logical in list.vertex.attributes.active.\n")
  }
  if((min(v,Inf) < 1) || (max(v,-Inf) > network.size(x))){  # combinatorial breakage:
    stop("Illegal vertex id specified in v in activate.edge.attribute.\n")
  }
    rule<-match.arg(rule)
    attributes <- list.vertex.attributes(x)
    attributes.active <- grep(".active",attributes,value=TRUE) 
    currently.active <- sapply(gsub(".active","",attributes.active),function(attr){
      any(!is.na(get.vertex.attribute.active(x,attr,onset=onset, terminus=terminus,length=length, at=at, rule = rule, na.omit = na.omit, active.default = active.default, dynamic.only=dynamic.only, require.active=require.active, return.tea=TRUE,null.na=TRUE)))
      })         
  if(dynamic.only){
   if (length(currently.active)>0){
     return(attributes.active[currently.active])
   } else {
     return(character(0))
   }
  } else {
    if (length(currently.active)>0){
      return(c(setdiff(attributes,attributes.active),attributes.active[currently.active]))
    } else {
      return(attributes)
    }
  }
  
}



##########################################################################
#list.edge.attributes.active
##########################################################################
list.edge.attributes.active <-function(x, onset=NULL, terminus=NULL,length=NULL, at=NULL, na.omit = FALSE, rule = c("any", "all"), e=seq_along(x$mel), require.active= FALSE, active.default = TRUE, dynamic.only=FALSE){
  # validate inputs
  
  if (!is.network(x)){
    stop("The first argument for get.edge.attributes.active must be a network object")
  }
  
  if(!is.null(at)) {
    if(!is.vector(at) || !is.numeric(at))
      stop("Singular time points given by the 'at' argument must be a numeric vector in list.edge.attributes.active\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be a numeric vector in list.edge.attributes.active\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric vector in list.edge.attributes.active\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric vector in list.edge.attributes.active\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(is.null(terminus) || is.null(length))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  rule<-match.arg(rule)
  # work around bug in network for 0 edges
  if (network.edgecount(x)>0)
  {
    attributes <- list.edge.attributes(x) 
    attributes.active <- grep(".active",attributes,value=TRUE)
    currently.active <- sapply(gsub(".active","",attributes.active),function(attr){
      any(!is.na(get.edge.value.active(x,attr,onset=onset, terminus=terminus,length=length, at=at, rule = rule, active.default = active.default, dynamic.only=dynamic.only, require.active=require.active,return.tea=TRUE)))
    }) 
  if(dynamic.only){
    if(length(currently.active)>0){
      return(attributes.active[currently.active])
    } else {
      return(character(0))
    }
  } else {
    if(length(currently.active)>0){
      return(c(setdiff(attributes,attributes.active),attributes.active[currently.active]))
    } else {
      return(attributes)
    }
  }
  } else { 
    # no edges to test for attributes
    return(character(0))
  }    
}
  




##########################################################################
#list.network.attribute.active
##########################################################################

list.network.attributes.active <-function(x, onset=NULL, terminus=NULL,length=NULL, at=NULL, na.omit = FALSE, rule = c("any", "all"), dynamic.only=FALSE){
  
  # checks for proper inputs
  if(!is.network(x)) 
    stop("list.network.attribute.active requires an argument of class network.\n")
  if(!is.null(at)) {
    if( !is.numeric(at))
      stop("'at' argument must be  numeric in list.network.attribute.active\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    if(length(at)!=1){
      stop("'at' time argument must have length one in list.network.attribute.active")
    }
  } else {
    if(!is.null(onset) &&  !is.numeric(onset))
      stop("Onset time must be  numeric  in list.network.attribute.active\n")
    if(!is.null(onset) &&  length(onset)!=1){
      stop("Onset time argument must have length one list.network.attribute.active\n")
    }
    if(!is.null(terminus) && !is.numeric(terminus))
      stop("Terminus times must be a numeric vector list.network.attribute.active\n")
    if(!is.null(terminus) &&  length(terminus)!=1){
      stop("Terminus time argument must have length one list.network.attribute.active\n")
    }
    if(!is.null(length) && (!is.numeric(length) || any(length < 0)))
      stop("Interval lengths must be a non-negative numeric in list.network.attribute.active\n")
    if(!is.null(length) &&  length(length)!=1){
      stop("Length time argument must have length one list.network.attribute.active\n")
    }
    
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(is.null(terminus) || is.null(length))
        stop("Query intervals must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
    rule<-match.arg(rule)
    attributes <- list.network.attributes(x)
    active.attributes <-grep(".active",attributes,value=T)
    currently.active <- sapply(gsub(".active","",active.attributes),function(attr){
      any(!is.na(get.network.attribute.active(x,attr,onset=onset, terminus=terminus,length=length, at=at,rule = rule, dynamic.only=dynamic.only,return.tea=TRUE)))
    })                                        
 
  if(dynamic.only){
    if (length(currently.active)>0){
      return(active.attributes[currently.active])
    } else {
      return(character(0))
    }
  } else {
    if (length(currently.active)>0){
      return(c(setdiff(attributes,active.attributes),active.attributes[currently.active]))
    } else {
      return(attributes)
    }
  }
 
  
}

