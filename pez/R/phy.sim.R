#' Simulate phylogenies
#'
#' Simulate phylogenies under pure birth/death or as a function of
#' trait evolution
#'
#' \code{sim.bd.tree} simulates a pure birth/death speciation
#' model. There are two important things to note: (1) speciation is
#' randomised before extinction, and only one thing can happen to a
#' lineage per timestep. (2) This code works well for my purposes, but
#' absurd parameter values can cause the function to crash.
#'
#' \code{sim.bd.tr.tree} is an extension of \code{sim.bd.tree}, and
#' all its caveats apply to it. It additionally simulated the
#' evolution of a trait under Brownain motion
#' (\code{tr.walk}). Species' speciation/extinction rates change
#' depending on whether they have a trait value similar to other
#' species (\code{sp.tr}, \code{ext.tr}). When a speciation event
#' happens, the two daughters split evenly about the ancestor's trait
#' value, taking values half-way to whatever the nearest species'
#' value is. To be precise: \eqn{$p(speciate)_i = speciate_i + sp.tr
#' \times min(trait distance)$}{p(speciate) = speciate +
#' sp.tr*min.tr.dist}, \eqn{$p(extinct)_i = exinction_i + ext.tr
#' \times min(trait distance)$}{p(extinct) = extinction +
#' ext.tr*min.tr.dist}, where \eqn{$i$}{i} denotes each species.
#'
#' \code{edge2phylo} is an internal function for the
#' \code{\link{sim.phy}} and \code{\link{sim.meta}} function families,
#' which may be of use to you. Check those functions' code for
#' examples of use.
#'
#' These functions are closely related to \code{\link{sim.meta}}; the
#' latter are extensions that simulate meta-community structure at the
#' same time.
#' 
#' @param speciate probability each species will speciate in each
#' time-step (0-1)
#' @param extinction probability each species will go extinct in each
#' time-step (0-1)
#' @param time.steps number of time-steps for simulation
#' @return \code{\link[ape:phylo]{phylo}} object with random
#' tip.labels, and trait values if using \code{sim.br.tr.tree}.
#' @author Will Pearse
#' @seealso sim.meta scape
#' @examples
#' tree <- sim.bd.phy(0.1, 0, 10)
#' plot(tree)
#' @rdname sim.phy
#' @name sim.phy
#' @importFrom stats runif
#' @export
sim.bd.phy <- function(speciate=0.1, extinction=0.025, time.steps=20){
    #Setup
    edge <- matrix(c(1,2,1,3), byrow=TRUE, ncol=2)
    species <- c(TRUE, TRUE)
    edge.length <- c(1,1)
    extinct <- numeric()

    #Keep the next iteration's edges separate to avoid running them through too early
    next.edge <- edge
    next.edge.length <- edge.length
    next.species <- species
    new.edge <- 4

    #Loop over timesteps
    for(i in seq(time.steps)){
        #Loop over edges
        for(j in seq(nrow(edge))){
            #Need a breaker to stop things doing multiple things per timestep...
            breaker <- TRUE
            #We only care about species, not nodes
            if(species[j] == TRUE){
                #Speciate?
                if(breaker & runif(1) <= speciate){
                    next.edge <- rbind(next.edge, matrix(c(edge[j,2], new.edge, edge[j,2], new.edge+1), byrow=TRUE, nrow=2))
                    next.species[j] <- FALSE
                    next.species <- append(next.species, c(TRUE, TRUE))
                    next.edge.length <- append(next.edge.length, c(1, 1))
                    new.edge <- new.edge+2
                    breaker <- FALSE
                } else next.edge.length[j] <- next.edge.length[j] + 1
                #Extinction?
                if(breaker & runif(1) <= extinction){
                    next.species[j] <- FALSE
                    extinct <- append(extinct, edge[j,2])
                    breaker <- FALSE
                }
            }
        }
        edge <- next.edge
        edge.length <- next.edge.length
        species <- next.species
    }

    #Turn into ape::phylo and return
    return(edge2phylo(edge, species, extinct, edge.length))
}

#' @param tr.range vector of length two specifying boundaries for
#' trait values (see notes); initial two species will be at the 25th
#' and 75th percentiles of this space. See also \code{tr.wrap}
#' @param sp.tr speciation rate's interaction with the minimum
#' distance between a species and the species most similar to it (see
#' details)
#' @param ext.tr extinction rate's interaction with the minimum
#' distance between a species and the species most similar to it (see
#' details)
#' @param tr.walk at each time-step a species not undergoing
#' speciation or extinction has its trait value drawn from a
#' distribution centered at its current value and with a standard
#' deviation set by this value. I.e., this is the rate of the Brownian
#' motion trait evolution.
#' @param tr.wrap whether to force species' trait values to stay
#' within the boundary defined by \code{tr.range}; default TRUE.
#' @author Will Pearse
#' @rdname sim.phy
#' @importFrom stats runif rnorm quantile
#' @export
sim.bd.tr.phy <- function(speciate=0.1, extinction=0.025, time.steps=20, tr.range=c(0,1), sp.tr=2, ext.tr=1, tr.walk=0.2, tr.wrap=TRUE){
    #Setup
    edge <- matrix(c(1,2,1,3), byrow=TRUE, ncol=2)
    species <- c(TRUE, TRUE)
    edge.length <- c(1,1)
    extinct <- extinct.traits <- numeric()
    traits <- unname(quantile(tr.range, c(0.25,0.75)))
    .wrap <- function(x, range){
        if(x < range[1]) return(range[1])
        if(x > range[2]) return(range[2])
        return(x)
    }
    
    #Keep the next iteration's edges separate to avoid running them through too early
    next.edge <- edge
    next.edge.length <- edge.length
    next.species <- species
    new.edge <- 4
    next.traits <- traits
    
    #Loop over timesteps
    for(i in seq(time.steps)){
        #Get minimum trait distances (--> alter speciation and extinction rates)
        min.dist <- abs(outer(traits, traits, `-`))
        diag(min.dist) <- NA
        #Be careful; dead species are all NA and so can cause warnings
        min.dist <- apply(min.dist, 2, function(x) if(all(is.na(x))) NA else min(x, na.rm=TRUE))
        
        #Loop over edges
        for(j in seq(nrow(edge))){
            #Need a breaker to stop things doing multiple things per timestep...
            breaker <- TRUE
            #We only care about species, not nodes
            if(species[j] == TRUE){
                #Speciate?
                if(breaker & runif(1) <= (speciate + sp.tr * min.dist[j])){
                    next.edge <- rbind(next.edge, matrix(c(edge[j,2], new.edge, edge[j,2], new.edge+1), byrow=TRUE, nrow=2))
                    next.species[j] <- FALSE
                    next.species <- append(next.species, c(TRUE, TRUE))
                    next.traits <- append(next.traits, c(traits[j]-0.5*min.dist[j],traits[j]+0.5*min.dist[j]))
                    next.edge.length <- append(next.edge.length, c(1, 1))
                    new.edge <- new.edge+2
                    breaker <- FALSE
                } else next.edge.length[j] <- next.edge.length[j] + 1
                #Extinction?
                if(breaker & runif(1) <= (extinction + ext.tr * min.dist[j])){
                    next.species[j] <- FALSE
                    extinct <- append(extinct, edge[j,2])
                    extinct.traits <- append(extinct.traits, next.traits[j])
                    breaker <- FALSE
                    next.traits[j] <- NA
                }
                #Brownian motion on trait if nothing else
                if(breaker & !is.na(next.traits[j]))
                    next.traits[j] <- rnorm(1, next.traits[j], sd=tr.walk)
                if(tr.wrap & !is.na(next.traits[j]))
                    next.traits[j] <- .wrap(next.traits[j], tr.range)
            }
        }
        edge <- next.edge
        edge.length <- next.edge.length
        species <- next.species
        traits <- next.traits
    }
    traits[is.na(traits)] <- extinct.traits

    #Turn into ape::phylo and return
    return(edge2phylo(edge, species, extinct, edge.length, traits))
}

#' @param edge a two-column matrix where the first column is the start
#' node, the second the destination, as in
#' \code{\link[ape:phylo]{phylo}$edge}
#' @param s which of the rows in the edge matrix represent
#' extant species
#' @param e which of the tips in the edge matrix are extinct
#' (DEFAULT: empty vector, i.e., none)
#' @param el a vector to be used to give edge.length to the
#' phylogeny (default NA, i.e., none)
#' @param t if given (default NA), a vector to be used for traits
#' (\code{$traits} slot) in the phylogeny
#' @author Will Pearse
#' @rdname sim.phy
#' @export
edge2phylo <- function(edge, s, e=numeric(0), el=NA, t=NULL){
    spp.no <- sort(c(edge[s,2], e))
    spp.edges <- edge[,2] %in% spp.no
    to.change <- matrix(0, nrow=nrow(edge), ncol=ncol(edge))
    for(i in seq_along(spp.no))
        to.change[edge >= spp.no[i]] <- to.change[edge >= spp.no[i]] + 1
    edge <- edge - to.change + length(spp.no)
    edge[spp.edges,2] <- seq_along(spp.no)
    tree <- list(edge=edge, tip.label=paste("r", order(spp.no), sep="_"), edge.length=el, Nnode=length(unique(edge[,1])), traits=t)
    class(tree) <- "phylo"
    return(tree)
}
