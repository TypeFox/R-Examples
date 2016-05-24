#  File networkDynamic/R/extract.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2012,2013 the statnet development team
######################################################################
# Contents:
#
# "%t%.network"
# network.extract
# network.dynamic.check
# network.collapse
# "%k%"
######################################################################

#Operator form for network.extract
"%t%"<-function(x,at){
  network.extract(x=x,at=at)
}


#Function to take a temporal extract/cross-section of a dynamically extended
#network object. 
#retain.all.verticies added if you want to not remove inactive verticies so networks stay the same size for comparison
#TODO: rewrite this using get.inducedSubgraph?
old.network.extract<-function(x,onset=NULL,terminus=NULL,length=NULL, at=NULL,
                          rule=c("any","all"),active.default=TRUE,retain.all.vertices=FALSE,trim.spells=FALSE){
  rule<-match.arg(rule)
  # determine which nodes/edges are active
  # nodes activity is straight forward
  # edge activity depends on the activity of the edge, but
  # also the activity of it's in/out (tail/head) nodes
  
  activeV<-is.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
                     e=NULL,v=seq_len(x%n%"n"), rule=rule, active.default=active.default)
  activeE=logical(0)
  if(length(x$mel)){
    activeE<-is.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
                       e=valid.eids(x),v=NULL, rule=rule, active.default=active.default)
    nullE <- sapply(x$mel, is.null)
    inV = sapply(x$mel[!nullE], "[[", "inl")  # in nodes of edges
    outV = sapply(x$mel[!nullE], "[[", "outl")  # out nodes of edges
    activeTH = sapply(1:length(inV), function(x){activeV[inV[x]] && activeV[outV[x]]})  # are head and tail active?
  }
  if(retain.all.vertices){
    newVid <-seq.int(x%n%"n")
  } else {
    newVid <- cumsum(activeV)
  }
  
  # Create network
  n<-ifelse(retain.all.vertices, x%n%"n", sum(activeV))
  #if(n==0)
  #  return(list()) # from before we allowed nets size 0!
  net<-network.initialize(n)
  # Set network-level attributes
  net$gal<-as.list(x$gal)
  net%n%"n"<-n
  net%n%"mnext"<-1
  if(is.bipartite(net) && net%n%'bipartite' > 0){
    net%n%"bipartite"<-newVid[net%n%"bipartite"]
  }
  # Set vertex-level attributes
  if(n>0){
    if(retain.all.vertices){
      net$val<-as.list(x$val)
    } else {
      net$val<-as.list(x$val[activeV]) # safer to do this with attribute methods... 
    }
  }  
  # Add edges
  if(length(activeE)){
    activeETH = activeE & activeTH
    if(any(activeETH)){
      tail<-as.list(lapply(x$mel[!nullE][activeETH],function(z){newVid[z$outl]}))
      head<-as.list(lapply(x$mel[!nullE][activeETH],function(z){newVid[z$inl]}))
      atl<-as.list(lapply(x$mel[!nullE][activeETH],"[[","atl"))
      nam<-as.list(lapply(x$mel[!nullE][activeETH],function(z){names(z$atl)}))
      add.edges(net,tail=tail,head=head,names.eval=nam,vals.eval=atl)
    }
  }
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
  if (trim.spells){
    # delete extra spell data on vertices edges and attributes that would be outside the query time range.
    # not sure how this should handle censoring
    # https://statnet.csde.washington.edu/trac/ticket/216
    
    # need to handle the case of 'at' (onset=terminus) queries seperately to avoid returning
    # spell matrix with the 'null' spell Inf,Inf
    if (onset==terminus){
      # remove all spells and replace with query spell
      deactivate.vertices(net)
      activate.vertices(net,onset=onset,terminus=terminus)
      deactivate.edges(net)
      activate.edges(net,onset=onset,terminus=terminus)
    } else {
      
      # deactivate spells on either side of query spell
      # also have to avoid setting spells to -Inf,-Inf or Inf,Inf
      if(onset!=-Inf){
        deactivate.vertices(net,onset=-Inf,terminus=onset)
        deactivate.edges(net,onset=-Inf,terminus=onset)
      }
      if (terminus!=Inf){
        deactivate.vertices(net,onset=terminus,terminus=Inf)
        deactivate.edges(net,onset=terminus,terminus=Inf)
      }
    }
    # if retain.all was used, mark retained vertices as inactive
    if(retain.all.vertices){
      deactivate.vertices(net,v=which(!activeV))
    }
    
    # if there are edges, trim edge attributes to range
    if (network.edgecount(net)>0){
      active.edge.attrs <-gsub(".active","",grep(".active",list.edge.attributes(net),value=TRUE))
      if(length(active.edge.attrs)>0){
        for(attr in active.edge.attrs){
          if (onset==terminus){
            # todo: need to worry about missing edges / unset attributes getting value out of whack?
            values<-get.edge.value.active(net,prefix=attr,at=onset)
            # todo: might be faster to just delete from list?
            deactivate.edge.attribute(net,attr,onset=-Inf,terminus=Inf)
            activate.edge.attribute(net,attr,values,at=onset)
          } else {
            if (onset!=-Inf){
              deactivate.edge.attribute(net,attr,onset=-Inf,terminus=onset)
            }
            if (terminus!=Inf){
              deactivate.edge.attribute(net,attr,onset=terminus,terminus=Inf)
            }
          }
        }
      }
    }
    # trim vertex attributes to range
    active.vertex.attrs <-gsub(".active","",grep(".active",list.vertex.attributes(net),value=TRUE))
    if(length (active.vertex.attrs)>0){
      for(attr in active.vertex.attrs){
        if(onset==terminus){
          values<-get.vertex.attribute.active(net,prefix=attr,at=onset)
          deactivate.vertex.attribute(net,attr,onset=-Inf,terminus=Inf)
          activate.vertex.attribute(net,prefix=attr,values,at=onset)
        } else {
          if(onset!=-Inf){
            deactivate.vertex.attribute(net,attr,onset=-Inf,terminus=onset)
          }
          if(terminus!=Inf){
            deactivate.vertex.attribute(net,attr,onset=terminus,terminus=Inf)
          }
        }
      }
    }
    # trim network attributes to range
    active.net.attrs <-gsub(".active","",grep(".active",list.network.attributes(net),value=TRUE))
    if(length(active.net.attrs)>0){
      for(attr in active.net.attrs){
        if (onset==terminus){
          value <-get.network.attribute.active(net,prefix=attr,at=onset)
          deactivate.network.attribute(net,attr,onset=-Inf,terminus=Inf)
          activate.network.attribute(net,attr,value,at=onset)
        } else {
          if(onset!=-Inf){
            deactivate.network.attribute(net,attr,onset=-Inf,terminus=onset)
          }
          if(terminus!=Inf){
            deactivate.network.attribute(net,attr,onset=terminus,terminus=Inf)
          }
        }
      }
    }
  } # end tripm.spells bock
  # update net.obs.period censoring info
  net.obs.period<-net%n%'net.obs.period'
  if(!is.null(net.obs.period)){
    # truncate the observations to the onset and terminus value
    obs<-net.obs.period$observations
    # subset to just spells that intersect query period
    
    obs<-obs[sapply(obs,function(ob){spells.overlap(ob,c(onset,terminus))})]
    
    
    if (length(obs)>0){
      # modify the onset of the first and terminus of the last but don't expand
      if(onset>obs[[1]][1]){
        obs[[1]][1]<-onset
      }
      if(terminus<obs[[length(obs)]][2]){
        obs[[length(obs)]][2]<-terminus
      }
    } else {
      # create a new spell
      #obs<-list(c(onset,terminus))
      # or should we instead return a null spell?
      obs<-list(c(Inf,Inf))
    }
    net.obs.period$observations<-obs
    net%n%'net.obs.period'<-net.obs.period
  }
  set.nD.class(net)
  
}

network.extract<-function(x,onset=NULL,terminus=NULL,length=NULL, at=NULL,
                               rule=c("any","all"),active.default=TRUE,retain.all.vertices=FALSE,trim.spells=FALSE){
  rule<-match.arg(rule)
  
  # determine which nodes/edges are active
  # nodes activity is straight forward
  # edge activity depends on the activity of the edge, but
  # also the activity of it's in/out (tail/head) nodes
  
  activeV<-is.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
                     e=NULL,v=seq_len(x%n%"n"), rule=rule, active.default=active.default)
  activeE=logical(0)
  if(length(x$mel)){
    activeE<-is.active(x=x,onset=onset,terminus=terminus,length=length,at=at,
                       e=valid.eids(x),v=NULL, rule=rule, active.default=active.default)
    nullE <- sapply(x$mel, is.null)
    inV = sapply(x$mel[!nullE], "[[", "inl")  # in nodes of edges
    outV = sapply(x$mel[!nullE], "[[", "outl")  # out nodes of edges
    activeTH = sapply(1:length(inV), function(x){activeV[inV[x]] && activeV[outV[x]]})  # are head and tail active?
  }
  if(retain.all.vertices){
    newVid <-seq.int(x%n%"n")
  } else {
    newVid <- cumsum(activeV)
  }
  
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
  
  # Create network
  n<-ifelse(retain.all.vertices, x%n%"n", sum(activeV))
  #if(n==0)
  #  return(list()) # from before we allowed nets size 0!
  net<-network.initialize(n)
  # Set network-level attributes
  gal<-as.list(x$gal)
  gal[['n']]<-n
  gal[["mnext"]]<-1
  if(is.bipartite(x) && gal[["bipartite"]] > 0){
    gal[["bipartite"]]<-newVid[gal[["bipartite"]]]
  }
  # update net.obs.period censoring info
  if(!is.null(gal[['net.obs.period']])){
    # truncate the observations to the onset and terminus value
    obs<-gal[['net.obs.period']]$observations
    # subset to just spells that intersect query period
    
    obs<-obs[sapply(obs,function(ob){spells.overlap(ob,c(onset,terminus))})]
    
    
    if (length(obs)>0){
      # modify the onset of the first and terminus of the last but don't expand
      if(onset>obs[[1]][1]){
        obs[[1]][1]<-onset
      }
      if(terminus<obs[[length(obs)]][2]){
        obs[[length(obs)]][2]<-terminus
      }
    } else {
      # create a new spell
      #obs<-list(c(onset,terminus))
      # or should we instead return a null spell?
      obs<-list(c(Inf,Inf))
    }
    gal[['net.obs.period']]$observations<-obs
  } 
  # assign all the network attribute changes to the network in a single call
  set.network.attribute(net,names(gal),gal)
  
  # Set vertex-level attributes
  if(n>0){
    if(retain.all.vertices){
      net$val<-as.list(x$val)
    } else {
      net$val<-as.list(x$val[activeV]) # safer to do this with attribute methods... 
    }
  } 
  
  # Add edges
  if(length(activeE)){
    activeETH = activeE & activeTH
    if(any(activeETH)){
      tail<-as.list(lapply(x$mel[!nullE][activeETH],function(z){newVid[z$outl]}))
      head<-as.list(lapply(x$mel[!nullE][activeETH],function(z){newVid[z$inl]}))
      atl<-as.list(lapply(x$mel[!nullE][activeETH],"[[","atl"))
      nam<-as.list(lapply(x$mel[!nullE][activeETH],function(z){names(z$atl)}))
      add.edges(net,tail=tail,head=head,names.eval=nam,vals.eval=atl)
    }
  }
  
  if (trim.spells){
    # delete extra spell data on vertices edges and attributes that would be outside the query time range.
    # not sure how this should handle censoring
    # https://statnet.csde.washington.edu/trac/ticket/216
    
    # need to handle the case of 'at' (onset=terminus) queries seperately to avoid returning
    # spell matrix with the 'null' spell Inf,Inf
    if (onset==terminus){
      # remove all spells and replace with query spell
      #deactivate.vertices(net)
      #activate.vertices(net,onset=onset,terminus=terminus)
      splmat<-list(matrix(c(onset,terminus),ncol=2))
      set.vertex.attribute(net,'active',splmat)
      set.edge.attribute(net,'active',splmat)
      #deactivate.edges(net)
      #activate.edges(net,onset=onset,terminus=terminus)
    } else {
      
      # deactivate spells on either side of query spell
      # also have to avoid setting spells to -Inf,-Inf or Inf,Inf
      if(onset!=-Inf){
        deactivate.vertices(net,onset=-Inf,terminus=onset)
        deactivate.edges(net,onset=-Inf,terminus=onset)
      }
      if (terminus!=Inf){
        deactivate.vertices(net,onset=terminus,terminus=Inf)
        deactivate.edges(net,onset=terminus,terminus=Inf)
      }
    }
    # if retain.all was used, mark retained vertices as inactive
    if(retain.all.vertices){
      deactivate.vertices(net,v=which(!activeV))
    }
    
    # if there are edges, trim edge attributes to range
    if (network.edgecount(net)>0){
      active.edge.attrs <-gsub(".active","",grep(".active",list.edge.attributes(net),value=TRUE))
      if(length(active.edge.attrs)>0){
        for(attr in active.edge.attrs){
          if (onset==terminus){
            # todo: need to worry about missing edges / unset attributes getting value out of whack?
            values<-get.edge.value.active(net,prefix=attr,at=onset)
            # todo: might be faster to just delete from list?
            deactivate.edge.attribute(net,attr,onset=-Inf,terminus=Inf)
            activate.edge.attribute(net,attr,values,at=onset)
          } else {
            if (onset!=-Inf){
              deactivate.edge.attribute(net,attr,onset=-Inf,terminus=onset)
            }
            if (terminus!=Inf){
              deactivate.edge.attribute(net,attr,onset=terminus,terminus=Inf)
            }
          }
        }
      }
    }
    # trim vertex attributes to range
    active.vertex.attrs <-gsub(".active","",grep(".active",list.vertex.attributes(net),value=TRUE))
    if(length (active.vertex.attrs)>0){
      for(attr in active.vertex.attrs){
        if(onset==terminus){
          values<-get.vertex.attribute.active(net,prefix=attr,at=onset)
          deactivate.vertex.attribute(net,attr,onset=-Inf,terminus=Inf)
          activate.vertex.attribute(net,prefix=attr,values,at=onset)
        } else {
          if(onset!=-Inf){
           deactivate.vertex.attribute(net,attr,onset=-Inf,terminus=onset)
          }
          if(terminus!=Inf){
           deactivate.vertex.attribute(net,attr,onset=terminus,terminus=Inf)
          }
       }
      }
    }
    # trim network attributes to range
    active.net.attrs <-gsub(".active","",grep(".active",list.network.attributes(net),value=TRUE))
    if(length(active.net.attrs)>0){
      for(attr in active.net.attrs){
        if (onset==terminus){
          value <-get.network.attribute.active(net,prefix=attr,at=onset)
          deactivate.network.attribute(net,attr,onset=-Inf,terminus=Inf)
          activate.network.attribute(net,attr,value,at=onset)
        } else {
          if(onset!=-Inf){
           deactivate.network.attribute(net,attr,onset=-Inf,terminus=onset)
          }
          if(terminus!=Inf){
           deactivate.network.attribute(net,attr,onset=terminus,terminus=Inf)
          }
        }
      }
    }
  } # end tripm.spells bock
  
  set.nD.class(net)

}

# internal Check spell matrix function to veryfiy consistency of spell matricies
# returns a logical vector where first element says if dimensions are ok, second if spells are ok
chkspellmat<-function(z){
  goodDims=TRUE
  goodSpells=TRUE
  if(!is.null(z)) {
    if(length(dim(z))!=2) {
      goodDims=FALSE
    } else {
      if(dim(z)[2]!=2) {
        goodDims=FALSE
      } else {
        if(NROW(z)==1){ #Single spell - either equal, ordered or Inf,Inf
          if(all(z==Inf)||(z[1,2]>=z[1,1]))
            goodSpells=TRUE
          else
            goodSpells=FALSE
        }else{       #Multiple spells - equal, ordered, non-overlapping
          if(all(z[,2]>=z[,1])&&all(z[-1,1]-z[-NROW(z),2]>=0))
            goodSpells=TRUE
          else
            goodSpells=FALSE
        }
      }
    }
  }
  c(goodDims, goodSpells)
} # end spell mat function

# checks that a tea attribute has appropriate structure, returning true if ok
chkTeaAttrOk <- function(x){
  if(!is.list(x)){
    return(FALSE)
  }
  if (length(x)==1 && is.na(x[[1]])){
    return(TRUE) # its not bad, its just not there
  }
  if(length(x)!=2){
    return(FALSE)
  }
  if(!is.list(x[[1]])){
    return(FALSE)
  }
  if(!is.matrix(x[[2]])){
    return(FALSE)
  }
  if(length(x[[1]])!=nrow(x[[2]])){
    return(FALSE)
  }
  if(!all(chkspellmat(x[[2]]))){
    return(FALSE)
  }  
  
  return(TRUE)
}

#Function to check dynamic consistency of network objects.  Not terribly safe,
#long-term.
network.dynamic.check<-function(x,verbose=TRUE, complete=TRUE){
  if(!is.network(x))
    stop("network.dynamic.check requires an object of class network.\n")
  pass.complete=TRUE
  if (complete) {
    if(network.size(x)==0){
      vertok<-logical(0)
    } else {
      #Check to ensure that vertex activity matrices are legit
      vertok<-sapply(x$val,function(y){chkspellmat(y$active)})
      if(verbose && any(!vertok[1,])) {
        cat("The dimensionality of the spell matrices is incorrect for vertex/vertices ")
        cat(which(!vertok[1,]))
        cat(".\n")
      }
      if(verbose && any(!vertok[2,])) {
        cat("The ordering of the spell matrices is incorrect for vertex/vertices ")
        cat(which(!vertok[2,]))
        cat(".\n")
      }
      vertok <- apply(vertok, 2, function(x){x[1] & x[2]})
      if(any(!vertok)) pass.complete=FALSE
    }
  
    # check vertex attribute activity
    vrtAttrOk <- rep(TRUE,length=network.size(x))
    dynVrtAttrs <- grep(".active",list.vertex.attributes(x),value=TRUE)
    for (attr in dynVrtAttrs){
      badVrts <- which(!sapply(get.vertex.attribute(x,attr,unlist=FALSE),chkTeaAttrOk))
      if (length(badVrts)>0){
        vrtAttrOk[badVrts] <-FALSE;
        if (verbose){
          cat(paste("Dynamic vertex attribute '",attr,"' is malformed for vertex ids",paste(badVrts,collapse=" "),"\n"))
        }
      }
    }
    
    if(network.edgecount(x)>0){
      #Check to ensure that edge activity matrices are OK
      active <- lapply(lapply(x$mel, "[[", "atl"), "[[", "active")
      edgeok <- sapply(active, chkspellmat)
      if(verbose && any(!edgeok[1,])) {
        cat("The dimensionality of the spell matrices is incorrect for edge(s) ")
        cat(which(!edgeok[1,]))
        cat(".\n")
      }
      if(verbose && any(!edgeok[2,])) {
        cat("The ordering of the spell matrices is incorrect for edge(s) ")
        cat(which(!edgeok[2,]))
        cat(".\n")
      }
      edgeok = apply(edgeok, 2, function(x){x[1] & x[2]})
      if(any(!edgeok)) pass.complete=FALSE
      
      # check edge attribute activity
      edgeAttrOk <- rep(TRUE,length=network.edgecount(x))
      dynEdgeAttrs <- grep(".active",list.edge.attributes(x),value=TRUE)
      for (attr in dynEdgeAttrs){
        badE <- which(!sapply(get.edge.value(x,attr,unlist=FALSE),chkTeaAttrOk))
        if (length(badE)>0){
          edgeAttrOk[badE] <-FALSE;
          if (verbose){
            cat(paste("Dynamic edge attribute '",attr,"' is malformed for edge ids",paste(badE,collapse=" "),"\n"))
          }
        }
      }
    } else {
      edgeok <-logical(0)
      edgeAttrOk <- logical(0)
    }
    
    # check network attribute activity
    netAttrOk <- TRUE
    dynNetAttrs <- grep(".active",list.network.attributes(x),value=TRUE)
    for (attr in dynNetAttrs){
      badN <- chkTeaAttrOk(get.network.attribute(x,attr,unlist=FALSE))
      if (!badN){
        netAttrOk <-FALSE;
        if (verbose){
          cat(paste("Dynamic network attribute '",attr,"' is malformed\n"))
        }
      }
    }
    
    # check net.obs.period
    if('net.obs.period'%in%list.network.attributes(x)){
      #.check.net.obs.period(x%n%'net.obs.period')
      net.obsOK <-tryCatch(
        {.check.net.obs.period(x%n%'net.obs.period')
         TRUE},
        error=function(e){if(verbose){
          cat(paste("Problem with net.obs.period attribute:",conditionMessage(e),'\n'))
          }
          return(FALSE)
        }
      )
      # check for out-of-bound spells
      # this only checks the maximal extend of the range, not internal gaps
      obs.bounds<-range(unlist((x%n%'net.obs.period')$observations))
      net.bounds<-range(get.change.times(x))
      if(obs.bounds[1]>net.bounds[1] | obs.bounds[2]<net.bounds[2]){
        net.obsOK<-FALSE
        if (verbose){
          cat(paste("Network activity was detected outside the range indicated by the net.obs.period 'observations' component."))
        }
      }
      
    } else {
      net.obsOK<-NULL
    }
    
  } # end complete block

  if(pass.complete) {  
    dyadok<-rep(F,length(x$mel))
    for(i in seq_along(x$mel)) {
      if(is.null(x$mel[[i]])) {
        dyadok[i]<-TRUE
      } else {
        y<-x$mel[[i]]
        ep<-c(y$outl,y$inl)
        act<-y$atl$active
        if (is.null(act)) {
          dyadok[i] <- TRUE
        } else {
          #Verify endpoint activity
          flag<-TRUE
          for(j in seq_len(nrow(act)))          #Check each spell...
            if(any(act[j,]!=c(Inf,Inf)))        #Not a placeholder spell...
              for(k in seq_along(ep))           #...against each endpoint
                if(flag) {
                  active.try = try(flag <- is.active(x=x,onset=act[j,1],terminus=act[j,2],v=ep[k],
                                   rule="all"), silent=T)
                  if(class(active.try)=="try-error") {
                    flag <- FALSE
                    warn.text="This check encountered an error in the is.active function."
                    if(!complete) warn.text=paste(warn.text, "Re-run check with 'complete=T'")
                    warning(warn.text)
                  } else {
                    flag <- active.try
                  }
                }
          dyadok[i]<-flag
        }
      }
    }
    if(verbose && any(!dyadok)){
      cat("Edges were found active where the endpoints where not in edge(s) ")
      cat(which(!dyadok))
      cat(".\n")
    }
  }

  #Return the results
  if(complete && pass.complete) {
    list(vertex.checks=vertok, edge.checks=edgeok, dyad.checks=dyadok,vertex.tea.checks=vrtAttrOk,edge.tea.checks=edgeAttrOk,network.tea.checks=netAttrOk,net.obs.period.check=net.obsOK)
  } else if (complete) {
    list(vertex.checks=vertok, edge.checks=edgeok)
  } else {
    list(dyad.checks=dyadok)
  }
}

#Operator form for network.collapse
"%k%"<-function(dnet,at){
  network.collapse(dnet=dnet,at=at)
}

# execute a network crossection, and then an attribute crossection. 
old.network.collapse <- function(dnet,onset=NULL,terminus=NULL, at=NULL, length=NULL,rule=c("any","all","earliest","latest"),active.default=TRUE,retain.all.vertices=FALSE,rm.time.info=TRUE,...){
  # check args
  if(missing(dnet) || !is.networkDynamic(dnet)){
    stop("network.collapse requires that the first argument be a networkDynamic object")
  }
  
  
  if(!is.null(at)) {
    if(!is.numeric(at))
      stop("Activation times must be numeric in activate.network.attribute.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be numeric in network.collapse.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric value in network.collapse.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval length must be a non-negative numeric value in network.collapse.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }
  
  
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
    stop("Onset times must precede terminus times in network.collapse\n")
  }
  rule<-match.arg(rule)
  exRule<-'any'  # network extract doesn't support the 'earliest' and 'latest' rules
  if(rule=='all') exRule<-'all'
  
  # we do not need to trim spells (which may be expensive) if we won't be returning the activity count info anyway
  net <- network.extract(dnet,onset=onset,terminus=terminus,rule=exRule,active.default=active.default, trim.spells=!rm.time.info,retain.all.vertices=retain.all.vertices)
  # collapse network level attributes
  knownNAttrs<-list.network.attributes(net)
  activeAttrs <- knownNAttrs[grep(".active",knownNAttrs)]
  for(attr in activeAttrs){
    net<-set.network.attribute(net,sub(".active","",attr),get.network.attribute.active(net,sub(".active","",attr),onset=onset,terminus=terminus,rule=rule))
    net<-delete.network.attribute(net,attr)
  }
  
  if(network.size(net)>0){
    # only mess with edges if there are edges
    if (network.edgecount(net)>0){
      # collapse edge level attributes
      knownEdgeAttrs <-list.edge.attributes(net) # listing is expensive, avoid doing it multiple times
      activeEdgeAttrs <- knownEdgeAttrs[grep(".active",knownEdgeAttrs)]
      
      for (attr in activeEdgeAttrs){
        net<-set.edge.attribute(net,sub(".active","",attr),get.edge.value.active(net,sub(".active","",attr),onset=onset,terminus=terminus,rule=rule))
        net<-delete.edge.attribute(net,attr)
      }
      # handle edge activity.count and activity.duration
      if (!rm.time.info){
        # Todo: what  deleted edge cases?
        eteas <-get.edge.value(net,"active",unlist=FALSE)
        hasspls <- !sapply(eteas,is.null)
        # if there are no edges, don't bother updating
        if ('activity.count'%in%knownEdgeAttrs){
          warning('Edge attribute "activity.count" already exists and was not updated')
        } else {
          # assume existing edges with no activity are defined by single spell -Inf, Inf
          # if active.default=FALSE, edges with no activity will allready be removed
          edge.counts <-rep.int(1,length(eteas))
          edge.counts[hasspls] <-sapply(eteas[hasspls],nrow)
          set.edge.attribute(net,'activity.count',edge.counts)
        }
        if ('activity.duration'%in%knownEdgeAttrs){
          warning('Edge attribute "activity.duration" already exists and was not updated')
        } else {
          
          # assume existing edges with no activity are defined by single spell -Inf, Inf
          #  if active.default=FALSE, edges with no activity will allready have been removed)
          edge.durations<-rep.int(Inf,length(eteas))
          
          edge.durations[hasspls] <-sapply(eteas[hasspls],function(spls){sum(spls[,2]-spls[,1])},simplify=FALSE)
          set.edge.attribute(net,'activity.duration',edge.durations)
        }
      }
      delete.edge.attribute(net,'active')
    } # end edges block
    
    # collapse vertex level attributes
    knownVAttrs<-list.vertex.attributes(net)  # avoid listing multiple times
    activeNodeAttrs <-knownVAttrs[grep(".active",list.vertex.attributes(net))]
    for (attr in activeNodeAttrs){
      # have to decide if result should be unlisted so as not to mangle
      net<-set.vertex.attribute(net,sub(".active","",attr),get.vertex.attribute.active(net,sub(".active","",attr),onset=onset,terminus=terminus,rule=rule))
      net<-delete.vertex.attribute(net,attr)
    }
    
    # handle vertex activity.count and activity.duration
    if (!rm.time.info){
      vteas <- get.vertex.attribute(net,"active",unlist=FALSE)
      hasspls <-!is.na(vteas)
      if ('activity.count'%in%knownVAttrs){
        warning('Vertex attribute "activity.count" already exists and was not updated')
      } else {
        # TODO: probably need active default paramter to handle this
        # default to 'always active' assuming 1 fake spell
        vertex.counts <- rep.int(1,length(vteas))
        vertex.counts[hasspls] <-sapply(vteas[hasspls],nrow)
        set.vertex.attribute(net,'activity.count',vertex.counts)
      }
      if ('activity.duration'%in%knownVAttrs){
        warning('Edge attribute "activity.duration" already exists and was not updated')
      } else {
        # TODO: probably need active default paramter to handle this
        # default to 'always active' assuming so assume Inf duration 
        vertex.durations <-rep.int(Inf,length(vteas))
        vertex.durations[hasspls] <-sapply(vteas[hasspls],function(spls){sum(spls[,2,drop=FALSE]-spls[,1,drop=FALSE])},simplify=FALSE)
        set.vertex.attribute(net,'activity.duration',vertex.durations)
      }
    }
    delete.vertex.attribute(net,'active')
  } # end zero vertices block
  
  if (rm.time.info){
    # remove net.obs.period if it exists
    delete.network.attribute(net,'net.obs.period')
  }
  
  class(net)<-'network'
  return(net)
}

# execute a network crossection, and then an attribute crossection. 
network.collapse <- function(dnet,onset=NULL,terminus=NULL, at=NULL, length=NULL,rule=c("any","all","earliest","latest"),active.default=TRUE,retain.all.vertices=FALSE,rm.time.info=TRUE,...){
  # check args
  if(missing(dnet) || !is.networkDynamic(dnet)){
    stop("network.collapse requires that the first argument be a networkDynamic object")
  }
  

  if(!is.null(at)) {
    if(!is.numeric(at))
      stop("Activation times must be numeric in activate.network.attribute.\n")
    if(!(is.null(onset) && is.null(terminus) && is.null(length)))
      stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
  } else {
    if(!is.null(onset) && (!is.vector(onset) || !is.numeric(onset)))
      stop("Onset times must be numeric in network.collapse.\n")
    if(!is.null(terminus) && (!is.vector(terminus) || !is.numeric(terminus)))
      stop("Terminus times must be a numeric value in network.collapse.\n")
    if(!is.null(length) && (!is.vector(length) || !is.numeric(length) || any(length < 0)))
      stop("Interval length must be a non-negative numeric value in network.collapse.\n")
    if(!is.null(onset)) {
      if(!xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    } else {
      if(xor(is.null(terminus),is.null(length)))
        stop("Spells must be specified by exactly 1 of {at, onset+terminus, onset+length, length+terminus}")
    }
  }

  
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
    stop("Onset times must precede terminus times in network.collapse\n")
  }
  rule<-match.arg(rule)
  exRule<-'any'  # network extract doesn't support the 'earliest' and 'latest' rules
  if(rule=='all') exRule<-'all'
  
  # we do not need to trim spells (which may be expensive) if we won't be returning the activity count info anyway
  net <- network.extract(dnet,onset=onset,terminus=terminus,rule=exRule,active.default=active.default, trim.spells=!rm.time.info,retain.all.vertices=retain.all.vertices)
  #tracemem(net)
  # collapse network level attributes
  knownNAttrs<-list.network.attributes(net)
  activeAttrs <- knownNAttrs[grep(".active",knownNAttrs)]
  newNames <- sub(".active","",activeAttrs)
  newVals <- lapply(newNames, function(name){get.network.attribute.active(net,name,onset=onset,terminus=terminus,rule=rule,unlist=FALSE)})
  # check if we are setting multiple or single values, and adjust nesting as needed
  if(length(newNames)==1){
    newNames<-newNames[[1]]
    newVals<-newVals[[1]]
  }
  
  # set multiple values at once
  set.network.attribute(net,newNames,newVals)
  
  
  if (rm.time.info){
    # remove net.obs.period if it exists
    activeAttrs<-c(activeAttrs,'net.obs.period')
  }
  # delete multiple values at once
  delete.network.attribute(net,activeAttrs)
  
  
  
#     for(attr in activeAttrs){
#       net<-set.network.attribute(net,sub(".active","",attr),get.network.attribute.active(net,sub(".active","",attr),onset=onset,terminus=terminus,rule=rule))
#       net<-delete.network.attribute(net,attr)
#     }
  
  if(network.size(net)>0){
    # only mess with edges if there are edges
    if (network.edgecount(net)>0){
      # collapse edge level attributes
      knownEdgeAttrs <-list.edge.attributes(net) # listing is expensive, avoid doing it multiple times
      activeEdgeAttrs <- knownEdgeAttrs[grep(".active",knownEdgeAttrs)]
      
      for (attr in activeEdgeAttrs){
        net<-set.edge.attribute(net,sub(".active","",attr),get.edge.value.active(net,sub(".active","",attr),onset=onset,terminus=terminus,rule=rule))
        net<-delete.edge.attribute(net,attr)
      }
      # handle edge activity.count and activity.duration
      if (!rm.time.info){
        # Todo: what  deleted edge cases?
        eteas <-get.edge.value(net,"active",unlist=FALSE)
        hasspls <- !sapply(eteas,is.null)
        # if there are no edges, don't bother updating
        if ('activity.count'%in%knownEdgeAttrs){
          warning('Edge attribute "activity.count" already exists and was not updated')
        } else {
          # assume existing edges with no activity are defined by single spell -Inf, Inf
          # if active.default=FALSE, edges with no activity will allready be removed
          edge.counts <-rep.int(1,length(eteas))
          edge.counts[hasspls] <-sapply(eteas[hasspls],nrow)
          set.edge.attribute(net,'activity.count',edge.counts)
        }
        if ('activity.duration'%in%knownEdgeAttrs){
          warning('Edge attribute "activity.duration" already exists and was not updated')
        } else {
          
           # assume existing edges with no activity are defined by single spell -Inf, Inf
          #  if active.default=FALSE, edges with no activity will allready have been removed)
          edge.durations<-rep.int(Inf,length(eteas))
          
          edge.durations[hasspls] <-sapply(eteas[hasspls],function(spls){sum(spls[,2]-spls[,1])},simplify=FALSE)
          set.edge.attribute(net,'activity.duration',edge.durations)
        }
      }
      delete.edge.attribute(net,'active')
    } # end edges block
  
    # collapse vertex level attributes
    knownVAttrs<-list.vertex.attributes(net)  # avoid listing multiple times
    activeNodeAttrs <-knownVAttrs[grep(".active",list.vertex.attributes(net))]
    newNodeNames<-sub(".active","",activeNodeAttrs)
    if (length(newNodeNames)>0){
      newNodeVals<-lapply(newNodeNames,function(name){get.vertex.attribute.active(net,name,onset=onset,terminus=terminus,rule=rule,unlist=FALSE)})
      # have to avoid situation where we accidentally call singular version with plural nested values
      if (length(newNodeNames)==1){
        newNodeVals<-unlist(newNodeVals,recursive=FALSE)
      }
      set.vertex.attribute(net,attrname=newNodeNames,value=newNodeVals)
  #     for (attr in activeNodeAttrs){
  #       # have to decide if result should be unlisted so as not to mangle
  #       net<-set.vertex.attribute(net,sub(".active","",attr),get.vertex.attribute.active(net,sub(".active","",attr),onset=onset,terminus=terminus,rule=rule))
  #       net<-delete.vertex.attribute(net,attr)
  #     }
      # also remove the vertex activity attribute
    }
    
    
    # handle vertex activity.count and activity.duration
    if (!rm.time.info){
      vteas <- get.vertex.attribute(net,"active",unlist=FALSE)
      hasspls <-!is.na(vteas)
      if ('activity.count'%in%knownVAttrs){
        warning('Vertex attribute "activity.count" already exists and was not updated')
      } else {
        # TODO: probably need active default paramter to handle this
        # default to 'always active' assuming 1 fake spell
        vertex.counts <- rep.int(1,length(vteas))
        vertex.counts[hasspls] <-sapply(vteas[hasspls],nrow)
        set.vertex.attribute(net,'activity.count',vertex.counts)
      }
      if ('activity.duration'%in%knownVAttrs){
        warning('Edge attribute "activity.duration" already exists and was not updated')
      } else {
        # TODO: probably need active default paramter to handle this
        # default to 'always active' assuming so assume Inf duration 
        vertex.durations <-rep.int(Inf,length(vteas))
        vertex.durations[hasspls] <-sapply(vteas[hasspls],function(spls){sum(spls[,2,drop=FALSE]-spls[,1,drop=FALSE])},simplify=FALSE)
        set.vertex.attribute(net,'activity.duration',vertex.durations)
      }
    }
    activeNodeAttrs<-c(activeNodeAttrs,'active')
    delete.vertex.attribute(net,activeNodeAttrs)
  } # end zero vertices block
  
  # we don't know what additional classes it might have in the future, so just clobber networkDynamic part
  class(net)<-class(net)[class(net)!="networkDynamic"]
  return(net)
}


# a function to return multiple static network arguments corresponding to collapsed slices
get.networks <- function(dnet, start=NULL, end=NULL, time.increment=NULL, onsets=NULL, termini=NULL,...){
  
  # check args
  if(missing(dnet) || !is.networkDynamic(dnet)){
    stop("network.collapse requires that the first argument be a networkDynamic object")
  }
  
  

    if(!is.null(onsets) && (!is.vector(onsets) || !is.numeric(onsets)))
      stop("Onset times must be a numeric vector \n")
    if(!is.null(termini) && (!is.vector(termini) || !is.numeric(termini)))
      stop("Terminus times must be a numeric vector \n")
    if(!is.null(time.increment) && !is.numeric(time.increment))
      stop("time.increment must be a non-negative numeric value\n")
    
  
  
  # arguments are null, try to guess from net.obs.period
  net.obs.period<-dnet%n%'net.obs.period'
  if (!is.null(net.obs.period) & is.null(onsets) & is.null(termini)){
    if (is.null(start)){
      start<-min(unlist(net.obs.period$observations))
    }
    if (is.null(end)){
      end<-max(unlist(net.obs.period$observations))
    }
    if (is.null(time.increment) && is.numeric(net.obs.period$time.increment)){
      time.increment<-net.obs.period$time.increment
    }
  }
  
  # if time increment is still missing, default to 1
  if (is.null(time.increment)){
    time.increment<-1
  }
  
  
  
  # create the lists of onsets and termini to slice by
  # check that either onsets and termini OR start, end, increment specified but not both
  
  if(!is.null(start) & !is.null(end) & !is.null(time.increment)){
    # USING START & END
    if (!is.null(onsets) | !is.null(termini)){
      stop("onsets & termini cannot be specified with start & end arguments\n")
    }
    # watch out for infs
    if(is.infinite(start) | is.infinite(end)){
      stop("start and end values must be finite")
    }
    onsets<-seq(from=start,to=end-time.increment,by=time.increment)
    termini<-seq(from=start+time.increment,to=end,by=time.increment)
  } else if (!is.null(onsets) & !is.null(termini)){
    # USING ONSETS & TERMINI
    if (length(onsets)!=length(termini)){
      stop('onsets and termini must have the same number of elements')
    }
    if(any(onsets>termini)){
      stop("Onset times must precede terminus times\n")
    }
  } else {
    stop("Unable to infer appropriate onsets and termini values for extracting networks, start & end, or onsets and termini parameters must be specified")
  }
  
  net.list <- lapply(seq_along(onsets), function(x){
    network.collapse(dnet,onset=onsets[x],terminus=termini[x],...)
  } )
  
  return(net.list)  
}
