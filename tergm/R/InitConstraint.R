#  File R/InitConstraint.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
#============================================================================
# This file contains the following 12 functions for initializing empty
# constraint lists (each prependend with "InitConstraint")
#         <edges>                   <outdegreedist>
#         <degrees>=<nodedegrees>   <bd>
#         <degreesTetrad>           <indegrees>
#         <degreesHexad>            <outdegrees>
#         <degreedist>              <hamming>
#         <indegreedist>            <observed>
#============================================================================

##########################################################################################
# Each of the <InitConstraint.X> functions accepts an existing constraint list, 'conlist',
# and to this adds an empty constraint list for term X; if any arguments are passed besides
# 'conlist", execution will halt.
#
# --PARAMETERS--
#   conlist: a list, presumably of constraints for other terms
#
# --RETURNED--
#   conlist: updated to include the initialized empty constraint list for term X
#
##########################################################################################

InitConstraint.atleast<-function(conlist, lhs.nw, nw=NULL, ...){
  if(is.null(nw)) stop("Formation constraint ``atleast'' requires a baseline network.",call.=FALSE)
  if(network.naedgecount(nw)) stop("Baseline network passed to formation constraint ``atleast'' may not have missing dyads.")
  conlist$atleast<-list(nw=nw)

  conlist$atleast$free.dyads <- function(){
    standardize.network(!nw) # invert and standardize the network
  }

  conlist
}
#ergm.ConstraintImplications("atleast", c())

InitConstraint.atmost<-function(conlist, lhs.nw, nw=NULL, ...){
  if(is.null(nw)) stop("Dissolution constraint ``atmost'' requires a baseline network.",call.=FALSE)
  if(network.naedgecount(nw)) stop("Baseline network passed to dissolution constraint ``atmost'' may not have missing dyads.")
  conlist$atmost<-list(nw=nw)

  conlist$atmost$free.dyads <- function(){
    standardize.network(nw)
  }

  conlist
}
#ergm.ConstraintImplications("atmost", c())
