#  File R/impute.network.list.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
impute.network.list <- function(nwl, imputers=c(), nwl.prepend=list(), nwl.append=list()){
  # TODO: Make it possible to write one's own imputers. E.g., break the following out into impute.network.stop(), impute.network.next(), etc..
  IMPUTERS <- c("next", "previous", "majority", "0", "1")
  imputers <- IMPUTERS[pmatch(imputers,IMPUTERS)]
  if(any(is.na(imputers))) stop("Unknown imputation option: ",sQuote(imputer),". Impute options must be a character vector containing one or more of ", dQuote("next"), ", ", dQuote("previous"), ", ", dQuote("majority"), ", ", dQuote("0"), ", and/or ", dQuote("1"), ".")
  
  for(imputer in imputers){
    nwl <- switch(imputer,
                  previous = {
                    nwl <- c(nwl.prepend, nwl)
                    nwl.NA <- sapply(nwl, network.naedgecount)>0 # Update which networks have missing dyads.
                    
                    for(t in seq_along(nwl)[-1])
                      if(nwl.NA[t]){
                        # Workaround for a bug in network (Ticket #80 in Trac)
                        na.el <-as.edgelist(is.na(nwl[[t]]))
                        na.eids <- apply(na.el, 1, function(e) get.edgeIDs(nwl[[t]], e[1],e[2], na.omit=FALSE))
                        nwl[[t]] <- delete.edges(nwl[[t]], na.eids)
                        nwl[[t]][na.el] <- nwl[[t-1]][na.el]
                      }
                    # Remove the prepended networks.
                    nwl[length(nwl.prepend)+seq_len(length(nwl)-length(nwl.prepend))]
                  },
                  `next` = {
                    nwl <- c(nwl, nwl.append)
                    nwl.NA <- sapply(nwl, network.naedgecount)>0
                    
                    for(t in rev(seq_along(nwl)[-length(nwl)]))
                      if(nwl.NA[t]){
                        # Workaround for a bug in network (Ticket #80 in Trac)
                        na.el <-as.edgelist(is.na(nwl[[t]]))
                        na.eids <- apply(na.el, 1, function(e) get.edgeIDs(nwl[[t]], e[1],e[2], na.omit=FALSE))
                        nwl[[t]] <- delete.edges(nwl[[t]], na.eids)
                        nwl[[t]][na.el] <- nwl[[t+1]][na.el]
                      }
                    # Remove the appended networks.
                    nwl[seq_len(length(nwl)-length(nwl.append))]
                  },
                  majority = {
                    lapply(nwl, function(y){
                      impute <- sign(network.edgecount(y,na.omit=TRUE)/network.dyadcount(y,na.omit=TRUE)-0.5)
                      if(impute==0){# If exact tie, can't impute.
                        y
                      }else{
                        # Workaround for a bug in network (Ticket #80 in Trac)
                      na.el <-as.edgelist(is.na(y))
                      na.eids <- apply(na.el, 1, function(e) get.edgeIDs(y, e[1],e[2], na.omit=FALSE))
                        impute <- impute > 0
                        y <- delete.edges(y, na.eids)
                        y[na.el] <- impute
                        y
                      }
                    })
                  },
                  `0` = {
                    lapply(nwl, function(y){
                      # Workaround for a bug in network (Ticket #80 in Trac)
                      na.el <-as.edgelist(is.na(y))
                        na.eids <- apply(na.el, 1, function(e) get.edgeIDs(y, e[1],e[2], na.omit=FALSE))
                      y <- delete.edges(y, na.eids)
                    })
                  },
                  `1` = {
                    lapply(nwl, function(y){
                      # Workaround for a bug in network (Ticket #80 in Trac)
                      na.el <-as.edgelist(is.na(y))
                      na.eids <- apply(na.el, 1, function(e) get.edgeIDs(y, e[1],e[2], na.omit=FALSE))
                      y <- delete.edges(y, na.eids)
                      y[na.el] <- 1
                      y
                    })
                  }
                  )
  }

  nwl.NA <- sapply(nwl, network.naedgecount)>0
     
  if("previous" %in% imputers && nwl.NA[1]) warning("Imputation option `previous' cannot impute dyads of the first network in the series.")
  if("next" %in% imputers && nwl.NA[length(nwl.NA)]) warning("Imputation option `next' cannot impute dyads of the last network in the series.")
  if("majority" %in% imputers && any(nwl.NA)) warning("Imputation option `majority' encountered an exact tie.")
  
  nwl  
}
impute.network.list <- function(nwl, imputers=c(), nwl.prepend=list(), nwl.append=list()){
  # TODO: Make it possible to write one's own imputers. E.g., break the following out into impute.network.stop(), impute.network.next(), etc..
  IMPUTERS <- c("next", "previous", "majority", "0", "1")
  imputers <- IMPUTERS[pmatch(imputers,IMPUTERS)]
  if(any(is.na(imputers))) stop("Unknown imputation option: ",sQuote(imputer),". Impute options must be a character vector containing one or more of ", dQuote("next"), ", ", dQuote("previous"), ", ", dQuote("majority"), ", ", dQuote("0"), ", and/or ", dQuote("1"), ".")
  
  for(imputer in imputers){
    nwl <- switch(imputer,
                  previous = {
                    nwl <- c(nwl.prepend, nwl)
                    nwl.NA <- sapply(nwl, network.naedgecount)>0 # Update which networks have missing dyads.
                    
                    for(t in seq_along(nwl)[-1])
                      if(nwl.NA[t]){
                        # Workaround for a bug in network (Ticket #80 in Trac)
                        na.el <-as.edgelist(is.na(nwl[[t]]))
                        na.eids <- apply(na.el, 1, function(e) get.edgeIDs(nwl[[t]], e[1],e[2], na.omit=FALSE))
                        nwl[[t]] <- delete.edges(nwl[[t]], na.eids)
                        nwl[[t]][na.el] <- nwl[[t-1]][na.el]
                      }
                    # Remove the prepended networks.
                    nwl[length(nwl.prepend)+seq_len(length(nwl)-length(nwl.prepend))]
                  },
                  `next` = {
                    nwl <- c(nwl, nwl.append)
                    nwl.NA <- sapply(nwl, network.naedgecount)>0
                    
                    for(t in rev(seq_along(nwl)[-length(nwl)]))
                      if(nwl.NA[t]){
                        # Workaround for a bug in network (Ticket #80 in Trac)
                        na.el <-as.edgelist(is.na(nwl[[t]]))
                        na.eids <- apply(na.el, 1, function(e) get.edgeIDs(nwl[[t]], e[1],e[2], na.omit=FALSE))
                        nwl[[t]] <- delete.edges(nwl[[t]], na.eids)
                        nwl[[t]][na.el] <- nwl[[t+1]][na.el]
                      }
                    # Remove the appended networks.
                    nwl[seq_len(length(nwl)-length(nwl.append))]
                  },
                  majority = {
                    lapply(nwl, function(y){
                      impute <- sign(network.edgecount(y,na.omit=TRUE)/network.dyadcount(y,na.omit=TRUE)-0.5)
                      if(impute==0){# If exact tie, can't impute.
                        y
                      }else{
                        # Workaround for a bug in network (Ticket #80 in Trac)
                      na.el <-as.edgelist(is.na(y))
                      na.eids <- apply(na.el, 1, function(e) get.edgeIDs(y, e[1],e[2], na.omit=FALSE))
                        impute <- impute > 0
                        y <- delete.edges(y, na.eids)
                        y[na.el] <- impute
                        y
                      }
                    })
                  },
                  `0` = {
                    lapply(nwl, function(y){
                      # Workaround for a bug in network (Ticket #80 in Trac)
                      na.el <-as.edgelist(is.na(y))
                        na.eids <- apply(na.el, 1, function(e) get.edgeIDs(y, e[1],e[2], na.omit=FALSE))
                      y <- delete.edges(y, na.eids)
                    })
                  },
                  `1` = {
                    lapply(nwl, function(y){
                      # Workaround for a bug in network (Ticket #80 in Trac)
                      na.el <-as.edgelist(is.na(y))
                      na.eids <- apply(na.el, 1, function(e) get.edgeIDs(y, e[1],e[2], na.omit=FALSE))
                      y <- delete.edges(y, na.eids)
                      y[na.el] <- 1
                      y
                    })
                  }
                  )
  }

  nwl.NA <- sapply(nwl, network.naedgecount)>0
     
  if("previous" %in% imputers && nwl.NA[1]) warning("Imputation option `previous' cannot impute dyads of the first network in the series.")
  if("next" %in% imputers && nwl.NA[length(nwl.NA)]) warning("Imputation option `next' cannot impute dyads of the last network in the series.")
  if("majority" %in% imputers && any(nwl.NA)) warning("Imputation option `majority' encountered an exact tie.")
  
  nwl  
}
