#######################################################################
# rEMM - Extensible Markov Model (EMM) for Data Stream Clustering in R
# Copyrigth (C) 2011 Michael Hahsler
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


## Simple Markov Chain for EMM
## to make it fast we don't use generics but use a smc_prefix

## Creator
SimpleMC <- function(size = 10L) {
    size <- as.integer(size)

    new("SimpleMC", 
	    unused = size:1,
	    top = size,
	    counts = matrix(0, ncol=size, nrow=size),
	    initial_counts = structure(rep(0, size),
		    names=rep(NA, size)))
}

smc_size <- function(x) {
    length(x@unused) - x@top
}

## state names are stored as names of the initial counts!!!
smc_states <- function(x, ...) {
    names(x@initial_counts)[!is.na(names(x@initial_counts))] 
}
    
smc_names2index <- function(x, state) {
    match(state, names(x@initial_counts))
}

smc_expand <- function(x) {
    old_size <- length(x@initial_counts)
    
    new_x <- SimpleMC(old_size*2L)
    new_x@counts[1:old_size, 1:old_size] <- x@counts

    new_x@initial_counts[1:old_size] <- x@initial_counts
    names(new_x@initial_counts)[1:old_size] <- names(x@initial_counts)
    new_x@top <- old_size + x@top
    new_x@unused[(old_size+1):length(new_x@unused)] <- x@unused

    new_x
}

smc_compact <- function(x) {
    keep <- !is.na(names(x@initial_counts))
    x@counts <- x@counts[keep, keep, drop=FALSE]
    x@initial_counts <- x@initial_counts[keep]
    x@top <- 0L
    x@unused <- rep(NA_integer_, length(x@initial_counts))

    x
}


smc_addState <- function(x, state) {

    state <- unique(state)

    if(any(smc_containsState(x, state))) {
	cont <- smc_containsState(x, state)
	warning("State ", state[cont], " already exists! (not added)")
	state <- state[-cont]
    }

    ## expand if necessary (we double the size)
    ## Note: it would make sense to implements shrink
    while(x@top < length(state)) x <- smc_expand(x)

    pos <- x@unused[(x@top-length(state)+1L):x@top]
    x@unused[(x@top-length(state)+1L):x@top] <- NA
    x@top <- x@top-length(state)

    names(x@initial_counts)[pos] <- state
        
    ## we don't need to clean up counts. removeState does it!
    ## we don't need to initialize x@initial_counts[pos] <- 0

    x
}

smc_removeState <- function(x, state) {
    pos <- smc_names2index(x, state)

    if(any(is.na(pos))) {
	missing <- is.na(pos)
	warning("State ", state[missing], " does not exist!")
	pos <- pos[!missing]
	state <- state[!missing]
	if(length(pos)<1) return(x)
    }

    ## clean matrix
    x@counts[pos,]<-0
    x@counts[,pos]<-0


    old_top <- x@top
    
    x@initial_counts[pos] <-0
    names(x@initial_counts)[pos] <- NA
    
    x@top <- x@top + length(pos)
    x@unused[(old_top+1):x@top] <- pos
    x
}

smc_addTransition <- function(x, from, to, w = 1) {
    pos_f <- smc_names2index(x, from)
    pos_t <- smc_names2index(x, to)
    
    if(length(pos_f) != length(pos_t)) stop("Lengths do not match!")
    if(length(pos_f) != length(w)) w <- rep.int(w[1], length(pos_f))
    
    missing <- is.na(pos_f) | is.na(pos_t)
    if(any(missing)) {
	warning("Some states do not exist (edge not added)!")
    }

    for(i in 1:length(pos_f)) 
	x@counts[pos_f[i], pos_t[i]] <- x@counts[pos_f[i], pos_t[i]] + w[i]
    x
}

smc_setTransitions <- function(x, from, to, value = 0) {
  pos_f <- smc_names2index(x, from)
  pos_t <- smc_names2index(x, to)
  
  if(length(pos_f) != length(pos_t)) stop("Lengths do not match!")
  if(length(pos_f) != length(value)) value <- rep.int(value[1], length(pos_f))
  
  missing <- is.na(pos_f) | is.na(pos_t)
  if(any(missing)) {
    warning("Some states do not exist (edge not added)!")
  }
  
  for(i in 1:length(pos_f)) 
    x@counts[pos_f[i], pos_t[i]] <- value[i]
  x
}



smc_removeSelfTransition <- function(x) {
    diag(x@counts) <- 0
    x
}

smc_removeTransition <- function(x, from, to) {
    pos_f <- smc_names2index(x, from)
    pos_t <- smc_names2index(x, to)

    #for(i in 1:length(pos_f)) 
    #	x@counts[pos_f[i], pos_t[i]] <- 0
    
    x@counts <- .Call("clear_matrix", x@counts, pos_f, pos_t)

    x
}

smc_mergeStates <- function(x, state) {
    pos <- smc_names2index(x, state)
    
    if(any(is.na(pos))) {
	missing <- is.na(pos)
	warning("Some states do not exist (merge might fail)!")
	state <- state[!missing]
	pos <- pos[!missing]
    }

    if(length(pos)<2) {
	warning("Need at least 2 states to merge!")
	return(x)
    }
    
    to <- pos[1]
    
    ## initial counts
    x@initial_counts[to] <- sum(x@initial_counts[pos])

    ## update matrix
    x@counts[to,] <- colSums(x@counts[pos,])
    x@counts[,to] <- rowSums(x@counts[,pos])

    ## remove additional states
    x <- smc_removeState(x, state[-1])
    x
}

smc_fade <- function(x, f) {
    x@counts <- x@counts * f
    x@initial_counts <- x@initial_counts * f
    x
}

smc_countMatrix <- function(x) {
    used <- !is.na(names(x@initial_counts)) 
    if(length(used)<1) return(matrix(ncol=0, nrow=0))
    
    m <- x@counts[used, used, drop = FALSE]
    dimnames(m) <- list(
	    names(x@initial_counts)[used], 
	    names(x@initial_counts)[used])

    m
}

smc_initialCounts <- function(x) {
    x@initial_counts[!is.na(names(x@initial_counts))]
}

smc_containsState <- function(x, state) {
    !is.na(match(state, names(x@initial_counts)))
}

smc_as.igraph <- function(x) {
    graph.adjacency(smc_countMatrix(x), weighted=TRUE)
}

## convert to graph (needs package graph!) 
smc_as.graph <- function(x) {
    if(!.installed("graph")) stop ("Package graph needed! Please install from Bioconductor.")
    graph::graphAM(adjMat=smc_countMatrix(x), edgemode= "directed", 
	    values=list(weight=1))
}



