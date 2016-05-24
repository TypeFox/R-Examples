#  File R/stergm.utils.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
# Keeps a list of "named" graphic devices.
#
# Usage: get.dev(name)
#
# If a graphic device with a given name exists, switch to it. If not,
# tries to grab an "unnamed" device and gives it a name. If none are
# available, creates a new device, switches to it, and "remembers" the
# name.

get.dev <- local({
  devs <- list()
  function(name){    
    if(is.null(devs[[name]]) || dev.set(devs[[name]])!=devs[[name]]){
      # Try to find an "unnamed" device to take over.
      free.devs <- setdiff(dev.list(),unlist(devs))
      if(length(free.devs))
        dev.set(free.devs[1])
      else
        dev.new()
      
      devs[[name]] <<- dev.cur()
    }else dev.set(devs[[name]])
    return(devs[[name]])
  }
})

# A customized level plot that uses cyan for positive values, red for
# negative, white for close to 0, black for exactly 0, and ensures
# that the scales on both sides are the same.
.my.levelplot <- function(m,levels=80,...){
  bound <- max(na.omit(c(abs(m))))

  requireNamespace('lattice', quietly=TRUE)
  lattice::levelplot(m, at=unique(c(seq(-bound,-.Machine$double.eps,length.out=levels/2+1),
                 seq(.Machine$double.eps,+bound,length.out=levels/2+1))),
            col.regions=c(hsv(h=0,s=seq(1,2/levels,length.out=levels/2),v=1),rgb(0,0,0),
              hsv(h=.5,s=seq(2/levels,1,length.out=levels/2),v=1)),
            ...)
}


# A wrapper around network.extract
# extracts the network at the specified time point and and attaches
# a network attribute representing that time point
# as well as numeric vector named "lasttoggle" representing the age of every (off-diagonal, 
# non-bipartite crossing) possible dyad in the network.
# updates: lasttoggle is NULL when duration.dependent is FALSE
network.extract.with.lasttoggle <- function(nwd, at, duration.dependent){
	nw <- network.extract(nwd, at=at)
  # check if the appropriate pid is defined, and if not, add it
  if (is.null(nwd%n%'vertex.pid')){
	  nw %v% "tergm_pid" <- which(is.active(nwd, at=at, v=seq_len(network.size(nwd))))
  }
	if(duration.dependent==1){
		lttails <- lapply(nw$mel, "[[", "outl")
		ltheads <- lapply(nw$mel, "[[", "inl")
		ltlts <- lapply(lapply(lapply(nw$mel, "[[", "atl"), "[[", 
						"active"), function(x) suppressWarnings(max(x[x <= at])))
		
		ltm <- if (is.bipartite(nw)) 
					m <- matrix(-Inf, nw %n% "bipartite", network.size(nw) - 
									nw %n% "bipartite")
				else m <- matrix(-Inf, network.size(nw), network.size(nw))
		for (i in seq_along(ltlts)) if (ltlts[[i]] != -Inf) {
				e <- c(lttails[[i]], ltheads[[i]])
				if (!all(e)) 
					next
				if (!is.directed(nw)) 
					e <- c(min(e), max(e))
				if (is.bipartite(nw)) 
					e[2] <- e[2] - nw %n% "bipartite"
				m[e[1], e[2]] <- ltlts[[i]] - 1
			}
		m[m == -Inf] <- round(-.Machine$integer.max/2)
		lasttoggle <-to.lasttoggle.matrix(m, is.directed(nw), is.bipartite(nw))
	} 
	else {  # non-duration dependent model
		lasttoggle <- NULL
	}
	
	nw %n% "time" <- at
	nw %n% "lasttoggle" <- lasttoggle
	nw
}


to.networkDynamic.lasttoggle <- function(nw){
  nwd <- nw
  if(!is.null(nw %n% "lasttoggle")){

    lt.edges <- ergm.el.lasttoggle(nw)

    lt.edges <- lt.edges[lt.edges[,3]>round(-.Machine$integer.max/2),,drop=FALSE] 

    # The +1 after lt.edges[,3] is important: lasttoggle is shifted by -1 relative to
    # networkDynamic (at least for now).
    if(nrow(lt.edges)) nwd <- deactivate.edges(nwd, onset=-Inf, terminus=lt.edges[,3]+1, e=apply(lt.edges[,1:2,drop=FALSE],1,function(e) get.edgeIDs(nw, e[1], e[2])))
  }
  nwd<-delete.network.attribute(nwd, "time")
  nwd<-delete.network.attribute(nwd, "lasttoggle")
  class(nwd) <- c("networkDynamic","network")
  #attr(nwd,"end") <- nw %n% "time"
  
  nwd
}

networkDynamic.apply.changes <- function(nwd, changes){
  
  # if there are no changes, just return the existing network
  if(nrow(changes)==0){
    return(nwd)
  }
  
  ## Add edges that were never present in the initial network.
  extant.edges <- as.edgelist(nwd)
  changed.edges <- unique(changes[,c("tail","head"),drop=FALSE])
  new.edges <- changed.edges[!(paste(changed.edges[,1],changed.edges[,2]) %in% paste(extant.edges[,1],extant.edges[,2])),,drop=FALSE]
  nwd <- add.edges(nwd,as.list(new.edges[,1]),as.list(new.edges[,2]))

  changes <- changes[order(changes[,"tail"], changes[,"head"], changes[,"time"], changes[,"to"]),,drop=FALSE]

  # Group changes by tail and head, resulting in a "data frame" with
  # columns tail, head, time, and to, with time and to being columns
  # of lists of changes for that dyad.
  changes <- aggregate(as.data.frame(changes[,c("time","to"),drop=FALSE]),by=list(tail=changes[,"tail"],head=changes[,"head"]),FUN=list)

  # Now, for each row in this "data frame", construct a list of lists,
  # each containing elements eID, times, and tos.
  changes <- apply(changes, 1, function(r)
                   list(eID=get.edgeIDs(nwd,r[["tail"]],r[["head"]]),
                        times=r[["time"]], tos=r[["to"]]))
  
  for(e in changes){
    tos <- e$tos
    times <- e$times
    eID <- e$eID
    
    if(!all(abs(diff(tos))==1)) stop("Problem with change matrix.")

    am <- nwd$mel[[eID]]$atl$active # Extant spells.
    
    if(tos[1]==0){ # First change is a dissolution.
      # No spell matrix:
      if(is.null(am)) am <- rbind(c(-Inf,+Inf))

      # If the last formation toggle is at the same time as the new
      # dissolution toggle, the whole spell gets dissolved away, with
      # the last row of am getting dropped below and not replaced by
      # anything. (This should not, normally, happen for the
      # simulate() functions.)
      #
      # Otherwise, prepend the onset of the extant tie.      
      if(am[nrow(am),1]==times[1]) times <- times[-1]
      else times <- c(am[nrow(am),1],times)

      # If ending with a formation, spell continues forever.
      if(tos[length(tos)]==1) times <- c(times, +Inf)

      # Construct a new spell matrix. (If times is empty, it's NULL, which still works.)
      am.new <- if(length(times)) matrix(times,ncol=2,byrow=TRUE)
      
      nwd[["mel"]][[eID]]$atl$active <- rbind(am[-nrow(am),],am.new)
    }else if(tos[1]==1){ # First change is a formation.

      # If the last dissolution toggle is at the same time as the new
      # formation toggle, the spell resumes as if never
      # dissolved. (This should not, normally, happen for the
      # simulate() functions.)
      if(!is.null(am) && am[nrow(am),2]==times[1]){
        times[1] <- am[nrow(am),1]
        am <- am[-nrow(am),,drop=FALSE]
      }

      # If ending with a formation, spell continues forever.
      if(tos[length(tos)]==1) times <- c(times, +Inf)
    
      # Construct a new spell matrix. (If times is empty, it's NULL, which still works.)
      am.new <- if(length(times)) matrix(times,ncol=2,byrow=TRUE)

      nwd[["mel"]][[eID]]$atl$active <- rbind(am,am.new)
    }
  }

  nwd
}
