#  File R/ergm.godfather.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
#=========================================================================
# This file contains the following 2 functions for computing changestat
# summaries of dynamic networks ??
#   <ergm.godfather>
#   <control.godfather>
#=========================================================================




###########################################################################
# <ergm.godfather>:  make the network a proposal it can't refuse. 
# Each toggle has a timestamp, and this function forces the network to make
# all of the changes at each unique timestamp value (in increasing order)
# keeping track of the change statistics that result. Thus, the final
# product is a matrix of change statistics in which the number of rows is
# determined by the # of unique timestamps and the number of columns is
# determined by the ERGM terms as usual.
#
# --PARAMETERS--
#   formula   : an ergm formula (i.e., nw ~ terms)
#   timestamps: a vector of timestamps for the given 'toggles'
#   toggles   : an edgelist of toggles that corresponds in length to
#               'timestamps'
#   sim       : a stergm sample, as returned by <stergm.getMCMCsample>
#               if 'sim' is not provided, both 'toggles' and
#               'timestamps' should be
#   start     : the start time; this is ignored if 'sim' is provided;
#               default=min(timestamps)
#   end       : the end time; this is ignored if 'sim' is provided;
#               default=max(timestamps)
#   accumulate: whether to proceed to the next timestamp without making
#               the proposed toggles (T or F); FALSE will force all toggles
#               to be realized on the network given in 'formula'
#               ?? if this is TRUE
#   verbose   : whether this and the C function should be verbose (T or F)
#               default=FALSE
#   control   : a list of additional tuning parameters for this function,
#               as returned by <control.godfather>;
#               default=<control.godfather>()
#
# --RETURNED--
#   the dynamic changestats summary as a list of the following:
#    stats     : a matrix, where the i,j entry represents the change in the
#                jth summary statistic between the original network and the
#                network at the ith unique timestamp
#    timestamps: the vector  c(NA, start:end)), where start and end are
#                specified either by attributes of 'sim' or by the 'start'
#                and 'end' inputs or default to the minimum and maximum
#                timestamps
#    newnetwork: the network after all toggles have been made if requested
#                by 'return_new_network' in <control.godfather>;
#                NULL otherwise
#
############################################################################

tergm.godfather <- function(formula, changes=NULL, toggles=changes[,-4,drop=FALSE],
                           start=NULL, end=NULL,
                           end.network=FALSE,
                           stats.start=FALSE,
                           verbose=FALSE,
                           control=control.tergm.godfather()){
  check.control.class("tergm.godfather")

  nw <- ergm.getnetwork(formula)
  
  formula <- ergm.update.formula(formula, nw~., from.new="nw")
  
  if(is.networkDynamic(nw)){
    if(!is.null(toggles)) stop("Network passed already contains change or toggle information.")

    toggles <- do.call(rbind, lapply(nw$mel, function(e) if(length(c(e$atl$active)[is.finite(c(e$atl$active))])) cbind(c(e$atl$active)[is.finite(c(e$atl$active))], e$outl,e$inl) else NULL))
    toggles[,1] <- ceiling(toggles[,1]) # Fractional times take effect at the end of the time step.
    
    net.obs.period<-nw%n%'net.obs.period'
    nwend <- if(!is.null(net.obs.period)) .get.last.obs.time(nw) else NULL
    nwstart<- if(!is.null(net.obs.period)) .get.first.obs.time(nw) else NULL
   
    start <- NVL(start,
                 nwstart,
                 suppressWarnings(min(toggles[,1]))-1
                 )
    if(start==Inf) stop("networkDynamic passed contains no change events or attributes. You must specify start explicitly.")

    end <- NVL(end,
               nwend,
               suppressWarnings(max(toggles[,1]))
               )
    if(end==-Inf) stop("networkDynamic passed contains no change events or attributes. You must specify end explicitly.")

    # The reason why it's > start is that the toggles that took effect
    # at start have already been applied to the network. Conversely,
    # it's <= end because we do "observe" the network at end, so we
    # need to apply the toggles that take effect then.
    toggles <- toggles[toggles[,1]>start & toggles[,1]<=end,,drop=FALSE]

    # This is important, since end is inclusive, but terminus is exclusive.
    if(!all(is.active(nw, onset=start, terminus=end+.Machine$double.eps*end*2, v=seq_len(network.size(nw)), rule="any")
            ==is.active(nw, onset=start, terminus=end+.Machine$double.eps*end*2, v=seq_len(network.size(nw)), rule="all")))
      stop("Network size and/or composition appears to change in the interval between start and end. This is not supported by ergm.godfather() at this time.")

    # Finally, we are ready to extract the network.
  duration.dependent <- if(is.durational(formula)){1} else {0}
  nw <- network.extract.with.lasttoggle(nw, at=start, duration.dependent)

  }else{
    if(is.null(toggles)) stop("Either pass a networkDynamic, or provide change or toggle information.")
      
    start <- NVL(start,
                 attr(toggles, "start"),
                 min(toggles[,1])-1
                 )
    end <- NVL(end,
               attr(toggles, "end"),
               max(toggles[,1])
               )

    # The reason why it's > start is that the toggles that took effect
    # at start have already been applied to the network. Conversely,
    # it's <= end because we do "observe" the network at end, so we
    # need to apply the toggles that take effect then.
    toggles <- toggles[toggles[,1]>start & toggles[,1]<=end,,drop=FALSE]
    
    if(is.null(nw %n% "lasttoggle")) nw %n% "lasttoggle" <- rep(round(-.Machine$integer.max/2), network.dyadcount(nw))
    nw %n% "time" <- start
  }

  if(!is.directed(nw)) toggles[,2:3] <- t(apply(toggles[,2:3,drop=FALSE], 1, sort))
  toggles <- toggles[order(toggles[,1],toggles[,2],toggles[,3]),,drop=FALSE]

  formula <- ergm.update.formula(formula, nw~., from.new="nw")
  m <- ergm.getmodel(formula, nw, expanded=TRUE, role="target")
  Clist <- ergm.Cprepare(nw, m)
  m$obs <- summary.statistics.network(m$formula)
  if(end.network){
    maxedges.sd <- sqrt(nrow(toggles)*0.25)*2 # I.e., if each toggle has probability 1/2 of being in a particular direction, this is the s.d. of the number of edges added.
    maxedges <- Clist$nedges + maxedges.sd*control$GF.init.maxedges.mul
  }

  if(verbose) cat("Applying changes...\n")
  repeat{
    z <- .C("godfather_wrapper",
            as.integer(Clist$tails), as.integer(Clist$heads),
            time = if(is.null(Clist$time)) as.integer(0) else as.integer(Clist$time),
            lasttoggle = as.integer(NVL(Clist$lasttoggle,0)),             
            as.integer(Clist$nedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.double(Clist$inputs),
            as.integer(nrow(toggles)), as.integer(toggles[,1]),
            as.integer(toggles[,2]), as.integer(toggles[,3]),
            as.integer(start), as.integer(end),
            s = double((1+end-start) * Clist$nstats),
            if(end.network) as.integer(maxedges) else as.integer(0),
            newnwtails = if(end.network) integer(maxedges+1) else integer(0),
            newnwheads = if(end.network) integer(maxedges+1) else integer(0),
            as.integer(verbose),
            status = integer(1), # 0 = OK, MCMCDyn_TOO_MANY_EDGES = 1
            PACKAGE="tergm")

    if(z$status==0) break;
    if(z$status==1){
      maxedges <- 5*maxedges
      if(verbose>0) message("Too many edges encountered in the simulation. Increasing capacity to ", maxedges)
    }
  }

  stats <- matrix(z$s + m$obs, ncol=Clist$nstats, byrow=TRUE)
  colnames(stats) <- m$coef.names
  if(!stats.start) stats <- stats[-1,,drop=FALSE]
  stats <- mcmc(stats, start=if(stats.start) start else start+1)
  
  if(end.network){ 
    if(verbose) cat("Creating new network...\n")
    newnetwork <- newnw.extract(nw,z)
    newnetwork %n% "time" <- z$time
    newnetwork %n% "lasttoggle" <- z$lasttoggle

    attr(newnetwork,"stats")<-stats
    newnetwork
  }else stats
}




####################################################################
# The <control.godfather> function allows for tuning of the
# <ergm.godfather> function
#
# --PARAMETERS--
#   maxedges          : the maximum number of edges to make space
#                       for for the new network; this is ignored
#                       if 5*Clist$nedges is greater; this is also
#                       ignored if 'return_new_network' is FALSE;
#                       default=100000
#
#
# --RETURNED--
#   a list of the above parameters
#
####################################################################

control.tergm.godfather<-function(GF.init.maxedges.mul=5
              ){
    control<-list()
    for(arg in names(formals(sys.function())))
      control[arg]<-list(get(arg))

    control <- set.control.class()
    control
  }
