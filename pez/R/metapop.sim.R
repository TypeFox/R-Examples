#' Simulate a meta-community (and its phylogeny)
#'
#' \code{sim.meta.comm} simulates species moving through a
#' metacommunity. At each time-step each cell's next abundance for
#' each species is \code{env.quality} - \code{current.abundance} +
#' \code{stochastic}, and a species gets as many chances to migrate in
#' each time-step as it has cells (the same cell could migrate
#' multiple times). I use a Poisson for everything because I don't
#' want half-species (these are individuals), and keeping everything
#' in Poisson makes it easier to compare the relative rates of
#' everything.
#' 
#' \code{sim.meta.phy.comm} As above, but with a simulation of
#' phylogeny as well - there are no additional extinction parameters,
#' since extinction happens as a natural consequence of ecological
#' interactions.
#' 
#' @param size the length and width of the meta-community in grid
#' cells
#' @param n.spp number of species
#' @param timesteps number of time-steps (each discrete)
#' @param p.migrate probability that a group of species in each grid
#' cell will migrate to another grid cell each timestep (i.e., 10
#' cells occuped by a species --> 10*\code{p.migrate} chance of
#' migration)
#' @param env.lam \eqn{$\lambda$}{lambda} value for Poisson
#' distribution used to distribute environmental quality; essentially
#' the carrying capacity (for each species separately) for that cell
#' @param abund.lam \eqn{$\lambda$}{lambda} value for Poisson
#' distribution used to distribute initial abundances and abundance
#' after migration
#' @param stoch.lam \eqn{$\lambda$}{lambda} value for Poisson
#' distribution of noise added to the next step abundance
#' calculation. With equal chance, this is taken as either a positive
#' or a negative number (see details if you're confused as to why this
#' is Poisson!)
#' @note \code{\link{scape}} is a much more sophisticated simulation
#' of the biogeography, but requires you to supply a phylogeny. You
#' pays your money, you makes your choice.
#' @return For \code{sim.meta.comm} a list with a species-site matrix
#' as the first slot, and the environment as the second. Rownames of
#' the site-species are the List with the x and y co-ordinates of the
#' simulation grid pasted together; colnames are arbitrary species
#' names. \code{sim.meta.comm}, a \code{\link{comparative.comm}}
#' object (since we have now simulated a phylogeny), with the same
#' naming convention for the site names.  phylogeny.
#' @author Will Pearse
#' @rdname sim.meta
#' @name sim.meta
#' @seealso \code{\link{sim.phy}} \code{\link{scape}}
#' @importFrom stats rpois rbinom
#' @export
sim.meta.comm <- function(size=10, n.spp=8, timesteps=10, p.migrate=0.05, env.lam=10, abund.lam=5, stoch.lam=1){
    #Setup environment and abundances
    env <- matrix(rpois(size^2, env.lam), nrow=size, ncol=size)
    abundance <- array(rpois(size^2*n.spp, abund.lam), dim=c(size,size,n.spp))
    
    #Loop over for timesteps
    for(i in seq(timesteps)){
        #Loop over species
        for(j in seq(n.spp)){
            #Calculate new abundance in each cell
            present <- abundance[,,j]>0
            stoch <- rpois(1,stoch.lam) * sample(c(1,-1),1)
            abundance[,,j][present] <- env[present] - abundance[,,j][present] + stoch
            #Migration
            # - randomly choose number of migration events, then choose cells for it to happen in (for simulation ease)
            # - remember the species may have died out...
            present <- abundance[,,j]>0
            if(sum(present) > 0){
                n.migrate <- sum(rbinom(sum(present),1,p.migrate))
                cells <- which(present, arr.ind=TRUE)
                for(k in seq(n.migrate)){
                    index <- cells[sample(seq(nrow(cells)), 1),]
                    index[1] <- index[1] + sample(c(-1,0,1),1)
                    index[2] <- index[2] + sample(c(-1,0,1),1)
                    #Not a wrapped world; you can fall off the edge!
                    # - if they've migrated where there are already species, let's just add everything to that (for simplicity, more than anything :p)
                    if(all(index<=size))
                        abundance[index] <- env[index] - rpois(1,abund.lam) + abundance[index]
                }
            }
        }
        #Clean up negative species etc. (probably not necessary...)
        abundance[abundance < 0] <- 0
    }
    
    #Format output and return
    comm <- apply(abundance, 3, unlist)
    t <- dim(abundance)
    t <- expand.grid(seq_len(t[1]), seq_len(t[2]))
    rownames(comm) <- paste(t[,1], t[,2], sep=".")
    env <- data.frame(as.numeric(env), row.names=rownames(comm))
    return(list(comm=comm, environment=env))
}

#' \code{sim.meta.phy.comm} simulate a (sort of) meta-community
#' 
#' @param p.speciate probabilty that, at each timestep, a species will
#' speciate. A species can only speciate, migrate, or reproduce if it
#' has individuals!
#' @return \code{sim.meta.phy.comm} \code{\link{comparative.comm}}
#' object that describes the data; note that the rownames of the
#' community object refer to the \code{row.column} of the data in the
#' simulated grid assemblages.
#' @rdname sim.meta
#' @author Will Pearse
#' @importFrom stats rpois rbinom runif
#' @export
sim.meta.phy.comm <- function(size=10, n.spp=8, timesteps=10, p.migrate=0.3, env.lam=10, abund.lam=5, stoch.lam=1, p.speciate=0.05){
    #Setup environment and abundances
    env <- matrix(rpois(size^2, env.lam), nrow=size, ncol=size)
    abundance <- array(rpois(size^2*n.spp, abund.lam), dim=c(size,size,n.spp))
    #Setup phylogeny
    edge <- matrix(c(rep(1,n.spp),seq(from=2,length.out=n.spp)), ncol=2)
    phy.species <- rep(TRUE, n.spp)
    edge.length <- rep(1, n.spp)
    phy.abund.lookup <- data.frame(phy=seq(n.spp),abund=seq(n.spp))
    new.edge <- n.spp+2
    
    #Loop over for timesteps
    for(i in seq(timesteps)){
        to.add.abundance <- array(NA,dim=c(size,size,0))
        #Loop over species
        for(j in seq(dim(abundance)[3])){
            #Calculate new abundance in each cell
            present <- abundance[,,j]>0
            stoch <- rpois(1,stoch.lam) * sample(c(1,-1),1)
            abundance[,,j][present] <- env[present] - abundance[,,j][present] + stoch
            #Migration
            # - randomly choose number of migration events, then choose cells for it to happen in (for simulation ease)
            # - remember the species may have died out...
            present <- abundance[,,j]>0
            if(sum(present) > 0){
                n.migrate <- sum(rbinom(sum(present),1,p.migrate))
                cells <- which(present, arr.ind=TRUE)
                for(k in seq(n.migrate)){
                    index <- unname(cells[sample(seq(nrow(cells)), 1),])
                    index[1] <- index[1] + sample(c(-1,0,1),1)
                    index[2] <- index[2] + sample(c(-1,0,1),1)
                    #Not a wrapped world; you can fall off the edge!
                    # - if they've migrated where there are already species, let's just add everything to that (for simplicity, more than anything :p)
                    if(all(index<=size))
                        abundance[,,j][index[1],index[2]] <- env[index[1],index[2]] - rpois(1,abund.lam) + abundance[,,j][index[1],index[2]]
                }
            }
            #Speciation
            # - only one shot per species per timestep
            if(sum(present)>0){
                curr.edge <- phy.abund.lookup$phy[j]
                if(runif(1) <= p.speciate){
                    #Handle phylogeny
                    edge <- rbind(edge, matrix(c(edge[curr.edge,2], new.edge, edge[curr.edge,2], new.edge+1), byrow=TRUE, nrow=2))
                    phy.species[curr.edge] <- FALSE
                    phy.species <- append(phy.species, c(TRUE, TRUE))
                    edge.length <- append(edge.length, c(1, 1))
                    new.edge <- new.edge+2
                    #Handle community data
                    to.add.abundance <- array(c(to.add.abundance, rep(0,size^2)), dim=dim(to.add.abundance)+c(0,0,1))
                    cells <- which(present, arr.ind=TRUE)
                    index <- unname(cells[sample(seq(nrow(cells)), 1),])
                    to.add.abundance[,,dim(to.add.abundance)[3]][index[1],index[2]] <- ceiling(abundance[,,j][index[1],index[2]]/2)
                    abundance[,,j][index[1],index[2]] <- ceiling(abundance[,,j][index[1],index[2]] / 2)
                    #Handle lookup
                    phy.abund.lookup$phy[j] <- nrow(edge)-1
                    phy.abund.lookup <- rbind(phy.abund.lookup, c(nrow(edge),nrow(phy.abund.lookup)+1))
                } else edge.length[curr.edge] <- edge.length[curr.edge] + 1
            }
        }
        #Clean up
        if(dim(to.add.abundance)[3] > 0)
            abundance <- array(c(abundance, to.add.abundance), dim=c(size,size,dim(to.add.abundance)[3]+dim(abundance)[3]))
        abundance[abundance < 0] <- 0
    }
    
    #Produce ape::phylo
    species <- seq(nrow(edge)) %in% phy.abund.lookup$phy
    tree <- edge2phylo(edge, species, el=edge.length)
    phy.abund.lookup$phy[order(phy.abund.lookup$phy)] <- tree$tip.label

    #Format output and return
    comm <- apply(abundance, 3, unlist)
    colnames(comm) <- phy.abund.lookup[,1]
    t <- dim(abundance)
    t <- expand.grid(seq_len(t[1]), seq_len(t[2]))
    rownames(comm) <- paste(t[,1], t[,2], sep=".")
    env <- data.frame(as.numeric(env), row.names=rownames(comm))
    return(comparative.comm(tree, comm, env=env, force.root=1))
}
