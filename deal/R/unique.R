## unique.R
## Author          : Claus Dethlefsen
## Created On      : Tue Jan 15 17:06:23 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Thu Jul 24 10:23:42 2003
## Update Count    : 68
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Susanne Gammelgaard Bøttcher, Claus Dethlefsen
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################

unique.networkfamily <- function(x,incomparables=FALSE,equi=FALSE,timetrace=FALSE,epsilon=1e-12,...) {
  ## returns a nwf with redundant networks removed
  ## nwf must be sorted
  ## equi=T: just one representative for each equivalence class (note
  ## that an equivalence class here is defined as all networks with
  ## the same score).

    ## Algorithm:
    ## create vector of scores
    ## create unique vector of scores
    ## for each unique score (=equivalence-class)
    ##     find all networks with this score
    ##     for each of these networks
    ##     check if it is already in the list. If not, put it in.
    nwf <- x
    
  if (timetrace) {t1 <- proc.time();cat("[Unique ")}
  n <- length(nwf)

  tab <- rep(NA,n)
  for (i in 1:n)
    tab[i] <- nwf[[i]]$score

  utab <- unique(tab) 

  if (equi) {
      ens <- abs(diff(tab)) < epsilon
      idx <- (1:(n-1))[!ens]
      if (!ens[n-1]) idx <- c(idx,n)
      utab <- tab[idx]
  
    nwl <- list()
    for (i in 1:length(utab))
      nwl[[i]] <- nwf[[(1:n)[tab==utab[i]][1]]]
  }
  else { ## more work to do

    nwl <- list(nwf[[1]])
    ntab <- c(nwf[[1]]$score)
    
    for (i in 2:length(nwf)) {
      try <- nwf[[i]]
      same <- nwl[(1:length(nwl))[ntab==c(try$score)]]
      jump <- FALSE
      if (length(same)>0) {
        for (j in 1:length(same))
          if (nwequal( try, same[[j]]))
            {
              jump <- TRUE
              break
            }
        if (!jump) {
          nwl <- c(nwl,list(try))
          ntab <- c(ntab,c(try$score))
        }
      }
      else {
        nwl <- c(nwl,list(try))
        ntab <- c(ntab,c(try$score))
      }

    }
  } # else
  
  class(nwl) <- "networkfamily"

  if (timetrace) {
    t2 <- proc.time()
    cat((t2-t1)[1],"]\n")
  }
  
  nwl
  }

nwequal <- function(nw1,nw2) {
  ## check if nw1 and nw2 has same DAG
  ## Output: (T/F)
  stopifnot(nw1$n==nw2$n) ## must have the same number of nodes.
  n <- nw1$n

  for (node in 1:n) {
      p1 <- nw1$nodes[[node]]$parents
      p2 <- nw2$nodes[[node]]$parents

    if ( length(p1) !=
         length(p2) )
      return(FALSE)
    N <- length(p1)
    if ( N>0 ) 
      if (!all(sort(p1)==sort(p2))) return(FALSE)
    }
  return(TRUE)
}

elementin <- function(nw,nwl) {
  ## is the network nw in the list nwl?
  n <- length(nwl)
  tab <- rep(NA,n)
  for (i in 1:n)
    tab[i] <- nwl[[i]]$score
  same <- nwl[(1:length(nwl))[tab==c(nw$score)]]
  if (!length(same)>0) return(FALSE)
  for (i in 1:length(same))
    if (nwequal(nw,same[[i]])) return(TRUE)
  return(FALSE)
}


