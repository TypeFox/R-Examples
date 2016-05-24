#  File R/combine.networks.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
all.same <- function(x){
  if(length(x)==0) return(TRUE)
  v0 <- x[1]
  for(v in x[-1]) if(!identical(v0,v)) return(FALSE)
  return(TRUE)
}

# create a single block-diagonal network by combining multible networks

combine.networks <- function(nwl, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), blockname="NetworkID", detect.edgecov=TRUE, standardized=FALSE, require.last.attr=TRUE){
  if(any(sapply(nwl, is.bipartite))) .combine.networks.bipartite(nwl=nwl, ignore.nattr=ignore.nattr, ignore.vattr=ignore.vattr, ignore.eattr=ignore.eattr, blockname=blockname, detect.edgecov=detect.edgecov, standardized=standardized, require.last.attr=require.last.attr)
  else .combine.networks.unipartite(nwl=nwl, ignore.nattr=ignore.nattr, ignore.vattr=ignore.vattr, ignore.eattr=ignore.eattr, blockname=blockname, detect.edgecov=detect.edgecov, standardized=standardized, require.last.attr=require.last.attr)
}


.combine.networks.unipartite <- function(nwl, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), blockname="NetworkID", detect.edgecov=TRUE, standardized=FALSE, require.last.attr=TRUE){
  if(any(diff(sapply(nwl, is.directed)))) stop("All networks must have the same directedness.")

  if(!standardized) nwl <- lapply(nwl, standardize.network)
  ns <- sapply(nwl, network.size)
  blks <- c(0, cumsum(ns))

  out <- network.initialize(sum(ns), directed=is.directed(nwl[[1]]))


  # Concatenate network attributes. If you run into what looks like a covariate matrix, combine it correctly.
  
  for(a in setdiff(Reduce(intersect,lapply(nwl, list.network.attributes)),
                          ignore.nattr)){ # I.e., iterate through common attributes.
    vl <- lapply(nwl, get.network.attribute, a, unlist=FALSE)

    # Here, try to autodetect covariate matrices and combine them.
    if(detect.edgecov
       && all(sapply(vl, is.matrix))
       && all(sapply(vl, nrow)==ns)
       && all(sapply(vl, ncol)==ns)
       && all.same(sapply(vl, mode))){

      # A logical vector that extracts off-diagonal element of the ns*ns matrix.


      offdiags <- unlist(lapply(ns, function(n) c(diag(1,n)==0)))
      # It doesn't matter what the "filler" elements are, as long as
      # adding them doesn't add another category and it's not NA. So,
      # what this does is as follows: grab the off-diagonal elements
      # of each covariate matrix, concatenate them into one vector,
      # remove the NAs, and take 0 (if it's present) or the minimum
      # value. (0 as filler can be helpful for sparse matrix
      # representations.)
      dummyvals <- na.omit(unlist(lapply(vl, "c"))[offdiags])
      dummyval <- if(0 %in% dummyvals) 0 else min(dummyvals)
      m <- matrix(dummyval, sum(ns), sum(ns))
      mode(m) <- mode(vl[[1]])
      
      for(b in seq_along(vl)){
        inds <- blks[b]+seq_len(ns[b])
        m[inds, inds] <- vl[[b]]
      }

      vl <- m
    }
    
    out <- set.network.attribute(out, a, vl)
  }

  # Concatenate vertex attributes.
  
  for(a in setdiff(Reduce(intersect,lapply(nwl, list.vertex.attributes)),
                          ignore.vattr)){ # I.e., iterate through common attributes.
    out <- set.vertex.attribute(out, a,
                                do.call(c, lapply(nwl, get.vertex.attribute, a, unlist=FALSE))
                                )
  }

  # Add ties and attributes

  for(b in seq_along(nwl)){
    el <- rbind(as.edgelist(nwl[[b]]),as.edgelist(is.na(nwl[[b]])))
    eids <- apply(el, 1, function(e) get.edgeIDs(nwl[[b]], e[1], e[2], na.omit=FALSE))

    vals <- lapply(nwl[[b]]$mel,"[[","atl")[eids]
    names <- lapply(vals, names)

    out <- add.edges(out, el[,1]+blks[b], el[,2]+blks[b], names.eval=names, vals.eval=vals)
  }

  # Finally, add a vertex attribute specifying the blocks

  out <- set.vertex.attribute(out, blockname, rep(seq_along(ns),ns))

  out
}


.combine.networks.bipartite <- function(nwl, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), blockname="NetworkID", detect.edgecov=TRUE, standardized=FALSE, require.last.attr=TRUE){
  if(!all(sapply(nwl, is.bipartite))) stop("This function operates only on bipartite networks.")
  if(any(sapply(nwl, is.directed))) stop("Bipartite directed networks are not supported at this time.")

  if(!standardized) nwl <- lapply(nwl, standardize.network)
  ns <- sapply(nwl, network.size)
  es <- sapply(nwl, "%n%", "bipartite")
  eblks <- c(0, cumsum(es))
  bip <- eblks[length(eblks)]
  ablks <- cumsum(c(bip, ns-es))

  out <- network.initialize(sum(ns), directed=is.directed(nwl[[1]]), bipartite=bip)

  # Concatenate network attributes. If you run into what looks like a covariate matrix, combine it correctly.
  
  for(a in setdiff(Reduce(intersect,lapply(nwl, list.network.attributes)),
                          ignore.nattr)){ # I.e., iterate through common attributes.
    vl <- lapply(nwl, get.network.attribute, a, unlist=FALSE)

    # Here, try to autodetect covariate matrices and combine them.
    if(detect.edgecov
       && all(sapply(vl, is.matrix))
       && all(sapply(vl, nrow)==es)
       && all(sapply(vl, ncol)==ns-es)
       && all.same(sapply(vl, mode))){

      # It doesn't matter what the "filler" elements are, as long as
      # adding them doesn't add another category and it's not NA. So,
      # what this does is as follows: grab the elements of each
      # covariate matrix, concatenate them into one vector, remove the
      # NAs, and take 0 (if it's present) or the minimum value. (0 as
      # filler can be helpful for sparse matrix representations.)
      dummyvals <- na.omit(unlist(lapply(vl, "c")))
      dummyval <- if(0 %in% dummyvals) 0 else min(dummyvals)
      m <- matrix(dummyval, sum(es), sum(ns-es))
      mode(m) <- mode(vl[[1]])
      
      for(b in seq_along(vl)){
        einds <- eblks[b]+seq_len(es[b])
        ainds <- ablks[b]+seq_len(ns[b]-es[b])
        m[einds, ainds-sum(es)] <- vl[[b]]
      }

      vl <- m
    }
    
    out <- set.network.attribute(out, a, vl)
  }

  # Concatenate vertex attributes.
  
  for(a in setdiff(Reduce(intersect,lapply(nwl, list.vertex.attributes)),
                   ignore.vattr)){ # I.e., iterate through common attributes.
    vl <- lapply(nwl, get.vertex.attribute, a, unlist=FALSE)
    v <- vector(mode(vl[[1]]), sum(ns))

    for(b in seq_along(vl)){
        einds <- eblks[b]+seq_len(es[b])
        ainds <- ablks[b]+seq_len(ns[b]-es[b])
        v[einds] <- vl[[b]][seq_len(es[b])]
        v[ainds] <- vl[[b]][es[b]+seq_len(ns[b]-es[b])]
    }
    
    out <- set.vertex.attribute(out, a, v)
  }

  # Add ties and attributes

  for(b in seq_along(nwl)){
    el <- rbind(as.edgelist(nwl[[b]]),as.edgelist(is.na(nwl[[b]])))
    eids <- apply(el, 1, function(e) get.edgeIDs(nwl[[b]], e[1], e[2], na.omit=FALSE))

    vals <- lapply(nwl[[b]]$mel,"[[","atl")[eids]
    names <- lapply(vals, names)

    out <- add.edges(out, el[,1]+eblks[b], el[,2]-es[b]+ablks[b], names.eval=names, vals.eval=vals)
  }

  # Finally, add a vertex attribute specifying the blocks

  out <- set.vertex.attribute(out, blockname, rep(rep(seq_along(ns),2),c(es,ns-es)))

  out
}
