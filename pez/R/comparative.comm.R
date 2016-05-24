#' Creates a community comparative ecology object, the basis of all
#' functions in pez
#'
#' Basic checking of whether the input data match up is performed; you
#' need only supply \code{comm} and \code{phy}, nothing else is
#' mandatory. You can manipulate the internals of
#' \code{comparative.comm}, or use the wrappers inside \code{pez} to
#' keep everything in order. Examples of these features are given
#' below; they are described in detailed at \code{\link{cc.manip}}.
#' 
#' @param phy phylogeny (in \code{\link[ape:phylo]{phylo}} format) of
#' species
#' @param comm community \code{matrix} (as used in
#' \code{\link{vegan}}) with species as columns and rows as
#' communities. Must contain \code{rownames} and \code{colnames}; NAs
#' are not checked for but probably unwise.
#' @param traits \code{data.frame} of species traits, with
#' \code{rownames} matching \code{comm}. Saved in the \code{data} slot
#' of the resulting \code{comparative.comm} object for compatibility
#' with \code{\link[caper:comparative.data]{comparative.data}}.
#' @param env \code{data.frame} of environmental data with
#' \code{rownames} matching \code{comm}
#' @param warn whether to warn if species/sites are dropped when
#' creating object (default: TRUE)
#' @param force.root if \code{phy} is unrooted, a \code{root.edge} of
#' force.root will be added (default: -1, which means this will never
#' happen). Rarely needed, rarely advisable.
#' @note \code{comparative.comm} is compatible with
#' \code{\link[caper:comparative.data]{comparative.data}}; this means
#' that the slot for species' trait data is called \code{data}. I
#' appreciate this is somewhat unwieldy, but hopefully you agree it is
#' helpful in the long-term.
#' @return comparative.comm object
#' @author Will Pearse
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
#' #Subset on species, then sites
#' data <- data[1:5,]
#' data <- data[,1:5]
#' #Site and species can be manipulated
#' species(data)
#' sites(data)[1:3] <- c("lovely", "invert", "sites")
#' #Other data can be viewed
#' trait.names(data)
#' env.names(data)
#' #Get assemblage phylogenies of all sites
#' assemblage.phylogenies(data)
#' #Do some manual manipulation of your objects (NOTE: $data for traits)
#' data$data$new.trait <- sample(letters, nrow(data$comm), replace=TRUE)
#' @importFrom ape is.rooted cophenetic.phylo
#' @importFrom ade4 scalewt
#' @seealso \code{\link{plot.comparative.comm}} \code{\link{cc.manip}} \code{link[caper:comparative.data]{comparative.data}}
#' @export
comparative.comm <- function(phy, comm, traits=NULL, env=NULL, warn=TRUE, force.root=-1){
  #Assertions and argument handling
  #Phylogeny
  if(!inherits(phy, "phylo")) 
    stop("'", deparse(substitute(phy)), "' not of class 'phylo'")
  if(! is.rooted(phy)){
    if(force.root != -1){
      phy$root.edge <- force.root
    } else {
      stop("'", deparse(substitute(phy)), "' is not rooted or has a basal polytomy.")
    }
  }
  if(any(duplicated(phy$tip.label))) stop('Duplicate tip labels present in phylogeny')
  #Community matrix
  if(! is.matrix(comm)) stop("'", deparse(substitute(comm)), "' must be an object of class 'matrix'.")
  if(is.null(colnames(comm))) stop("'", deparse(substitute(comm)), "' must have column names (preferably species!)")
  if(is.null(rownames(comm))) stop("'", deparse(substitute(comm)), "' must have row names")
  if(any(is.na(comm))) stop("'", deparse(substitute(comm)), "' contains NAs")
  #Traits
  if(!is.null(traits)){
    if(! is.data.frame(traits)) stop("'", deparse(substitute(traits)), "' must be an object of class 'data.frame'.")
    if(is.null(colnames(traits))) stop("'", deparse(substitute(traits)), "' must have row names (preferably species!)")
    if(is.null(rownames(traits))) stop("'", deparse(substitute(traits)), "' must have row names")
  }
  #Environment
  if(!is.null(env)){
    if(! is.data.frame(env)) stop("'", deparse(substitute(env)), "' must be an object of class 'data.frame'.")
    if(is.null(colnames(env))) stop("'", deparse(substitute(env)), "' must have row names (preferably sites!)")
    if(is.null(rownames(env))) stop("'", deparse(substitute(env)), "' must have row names")
  }
  
  #Create intersection/drop lists and warn
  species.to.keep <- intersect(phy$tip.label, colnames(comm))
  if(!is.null(traits)) species.to.keep <- intersect(species.to.keep, rownames(traits))
  if(!is.null(env)) sites.to.keep <- intersect(rownames(comm), rownames(env)) else sites.to.keep <- rownames(comm)
  if(warn){
    if(length(setdiff(phy$tip.label, species.to.keep)) > 0)
      warning("Mismatch between phylogeny and other data, dropping ", length(setdiff(phy$tip.label, species.to.keep)), " tips")
    if(length(setdiff(colnames(comm), species.to.keep)) > 0)
      warning("Mismatch between community matrix and other data, dropping ", length(setdiff(colnames(comm), species.to.keep)), " columns")
    if(!is.null(traits) & length(setdiff(rownames(traits), species.to.keep)) > 0)
      warning("Mismatch between traits and other data, dropping ", length(setdiff(rownames(traits), species.to.keep)), " columns")
    if(length(setdiff(rownames(comm), sites.to.keep)) > 0)
      warning("Mismatch between community matrix and data, dropping ", length(setdiff(rownames(comm), sites.to.keep)), " rows")
    if(!is.null(env) & length(setdiff(rownames(env), sites.to.keep)) > 0)
      warning("Mismatch between env. data and other data, dropping ", length(setdiff(rownames(env), sites.to.keep)), " rows")
  }
  
  #Subset data and keep record
  # - make sure to re-order factor levels where necessary
  comm.sp.lost <- setdiff(colnames(comm), species.to.keep)
  comm <- comm[, colnames(comm) %in% species.to.keep]
  comm.sites.lost <- setdiff(rownames(comm), sites.to.keep)
  comm <- comm[rownames(comm) %in% sites.to.keep, ]
  if(any(dim(comm)==0))
    stop("ERROR: community data has no sites in common with rest of data")
  phy.sp.lost <- setdiff(phy$tip.label, species.to.keep)
  if(length(phy.sp.lost) == length(phy$tip.label))
    stop("ERROR: phylogeny has no species in common with rest of data")
  phy <- drop_tip(phy, phy.sp.lost)
  if(!is.null(traits)){
    traits.sp.lost <- setdiff(rownames(traits), species.to.keep)
    traits <- traits[rownames(traits) %in% species.to.keep, , drop = FALSE]
    if(any(dim(traits)==0))
      stop("ERROR: trait data has no species in common with rest of data")
  } else traits.sp.lost <- character(0)
  if(!is.null(env)){
    env.sites.lost <- setdiff(rownames(env), sites.to.keep)
    env <- env[rownames(env) %in% sites.to.keep, , drop = FALSE]
    if(any(dim(env)==0))
      stop("ERROR: environmental data has no sites in common with rest of data")
  } else env.sites.lost <- character(0)
  
  #Put species and sites in same order(s)
  # - leave phy alone (if caper needs it altered it will alter it)
  # - match everything to phylogeny's order (makes later work easier)
  comm <- comm[order(rownames(comm)), ]
  comm <- comm[, match(phy$tip.label, colnames(comm))]
  if(!is.null(traits)){
    traits <- traits[match(phy$tip.label, rownames(traits)), , drop = FALSE]
    traits <- traits[, colnames(traits), drop = FALSE]
  }
  if(!is.null(env)){
    env <- env[match(rownames(comm), rownames(env)), , drop = FALSE]
    env <- env[, colnames(env), drop = FALSE]
  }
  
  #Makee output, and return
  output <- list(phy=phy, comm=comm, data=traits, env=env, dropped=list(comm.sp.lost = comm.sp.lost,
                                                                 comm.sites.lost=comm.sites.lost,
                                                                 phy.sp.lost=phy.sp.lost,
                                                                 traits.sp.lost=traits.sp.lost,
                                                                 env.sites.lost=env.sites.lost), names=names)
  class(output) <- c("comparative.comm", "comparative.data")
  return(output)
}

#' @param x \code{comparative.comm} object to be printed
#' @param ... ignored
#' @rdname comparative.comm
#' @method print comparative.comm
#' @export
print.comparative.comm <- function(x, ...){
    #Argument checking
    if(!inherits(x, "comparative.comm"))
        stop("'", substitute(deparse(x)), "' not of class 'comparative.comm'")
    
    # basic summary data
    cat("Comparative community dataset of", ncol(x$comm), "taxa:\n")
    cat("Phylogeny:\n")
    cat("   ", x$phy$Nnode, " internal nodes\n", sep='')
    cat("Community data:\n")
    cat("   ", nrow(x$comm), " sites, ", ncol(x$comm), " taxa\n")
    
    if(!is.null(x$data)){
	    cat("Trait data:\n")
    	cat("   ", ncol(x$data), " variables\n")
    } else cat("Trait data: None\n")
    
    if(!is.null(x$env)){
    	cat("Environmental data:\n")
    	cat("   ", nrow(x$env), " sites, ", ncol(x$env), " variables\n")
    } else cat("Environmental data: None\n")
  
	# report on mismatch on merge
    if(length(x$dropped$phy.sp.lost)>0 | length(x$dropped$comm.sp.lost)>0 | length(x$dropped$data.sp.lost)>0){
  	 cat('Dropped taxa:\n')
  	 cat('   phy : ', length(x$dropped$phy.sp.lost), '\n')
  	 cat('   comm : ', length(x$dropped$comm.sp.lost), '\n')
     if(!is.null(x$data))
        cat('   trts : ', length(x$dropped$data.sp.lost), '\n')
    }
  if(length(x$dropped$comm.sites.lost)>0 | length(x$dropped$env.sites.lost)>0){
  	 cat('Dropped sites:\n')
       cat('   comm : ', length(x$dropped$comm.sites.lost), '\n')
  	 if(!is.null(x$env))
       cat('   env : ', length(x$dropped$env.sites.lost), '\n')
  }
}

#' Manipulating and examining comparative.comm objects
#'
#' As described in the vignette, we recommend using these wrappers to
#' manipulate species and site data, as it guarantees that everything
#' will be kept consistent across all parts of the
#' \code{\link{comparative.comm}} object. With them, you can drop
#' species, sites, and work directly with each part of your data. You
#' can also manipulate your \code{\link{comparative.comm}} object's
#' \code{phy}, \code{data}, \code{env}, and \code{comm} slots directly
#' if you wish, but altering the object directly yourself runs the
#' risk of things getting unsynchronised.
#' 
#' @param x \code{comparative.comm} object
#' @param sites numbers of sites to be kept or dropped from \code{x};
#' must be given as numbers. For example, \code{x[1:5,]}, or
#' \code{x[-1:-5,]}, but not \code{x[c("site a", "site b"),]}.
#' @param spp numbers of species to be kept or dropped from \code{x};
#' must be given as numbers. For example, \code{x[,1:5]}, or
#' \code{x[,-1:-5]}, but not \code{x[c("sp a", "sp b"),]}.
#' @param warn whether to warn if species/sites are dropped when
#' creating object (default: TRUE)
#' @param value when altering a \code{\link{comparative.comm}}
#' object's internal structure, the thing that you're inserting into
#' it!
#' @note As described in \code{\link{comparative.comm}}, each
#' \code{\link{comparative.comm}} object contains a phylogeny
#' (\code{$phy}) and a site-by-species community matrix (as used in
#' \code{\link{vegan}}). Optionally, it may contain a
#' \code{data.frame} of trait data (each row a species, each column a
#' trait ) *called \code{data}* for compatibility with
#' \code{\link[caper:comparative.data]{comparative.data}}.
#' @rdname cc.manip
#' @name cc.manip
#' @seealso comparative.comm plot.comaparative.comm
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
#' #Subset on species, then sites
#' data <- data[1:5,]
#' data <- data[,1:5]
#' #Site and species can be manipulated
#' species(data)
#' sites(data)[1:3] <- c("lovely", "invert", "sites")
#' #Other data can be viewed
#' trait.names(data)
#' env.names(data)
#' #Get assemblage phylogenies of all sites
#' assemblage.phylogenies(data)
#' #Add some trait/env data in
#' traits(data)$new.trait <- sample(letters, nrow(comm(data)), replace=TRUE)
#' env(data)$new.env <- sample(letters, ncol(comm(data)), replace=TRUE)
#' #Manipulate/check phylogeny and community matrix
#' phy(data) #...tree(data) works too...
#' comm(data)[1,3] <- 3
#' comm(data) <- comm(data)[-3,]
#' @export
"[.comparative.comm" <- function(x, sites, spp, warn=FALSE) {
  #Assertions and setup
  if(!inherits(x, "comparative.comm"))
    stop("'", substitute(deparse(x)), "' not of class 'comparative.comm'")    
  
  #Handle species
  if(!missing(spp)){
    if(is.null(spp)) stop("Null indices not permitted on comparative community data objects")
    if(is.numeric(spp) | is.logical(spp))
        spp.to.keep <- colnames(x$comm)[spp] else spp.to.keep <- intersect(species(x),spp)
    comm <- x$comm[, spp.to.keep]
    phy <- drop_tip(x$phy, setdiff(x$phy$tip.label, spp.to.keep))
    if(!is.null(x$data))
      traits <- x$data[spp.to.keep, ] else traits <- NULL
    new.x <- comparative.comm(phy, comm, traits, x$env, warn=warn)
  } else new.x <- x
  
  #Handle sites
  if(!missing(sites)){
    if(is.null(sites)) stop("Null indices not permitted on comparative community data objects")
    if(is.numeric(sites) | is.logical(sites))
        sites.to.keep <- rownames(x$comm)[sites] else sites.to.keep <- intersect(sites(x),sites)
    comm <- new.x$comm[sites.to.keep, ]
    if(!is.null(x$env))
      env <- new.x$env[sites.to.keep, ] else env <- new.x$env
    new.x <- comparative.comm(x$phy, comm, new.x$data, env, warn=FALSE)
  }
  
  #Warn of dropped species/sites (if asked)
  if(warn){
    orig.species <- colnames(x$comm)
    orig.sites <- rownames(x$comm)
    if(length(setdiff(orig.species, new.x$phy$tip.label)) > 0)
      warning("Mismatch between phylogeny and other data, dropping ", length(setdiff(orig.species, new.x$phy$tip.label)), " tips")
    if(length(setdiff(orig.sites, colnames(new.x$comm))) > 0)
      warning("Mismatch between community matrix and other data, dropping ", length(setdiff(orig.species, colnames(new.x$comm))), " columns")
    if(!is.null(new.x$data) & length(setdiff(orig.species, rownames(new.x$data))) > 0)
      warning("Mismatch between traits and other data, dropping ", length(setdiff(orig.species, rownames(new.x$data))), " columns")
    if(length(setdiff(rownames(new.x$comm), orig.sites)) > 0)
      warning("Mismatch between community matrix and data, dropping ", length(setdiff(orig.sites, rownames(new.x$comm))), " rows")
    if(!is.null(new.x$env) & length(setdiff(rownames(new.x$env), orig.sites)) > 0)
      warning("Mismatch between env. data and other data, dropping ", length(setdiff(orig.sites, rownames(new.x$env))), " rows")
  }
	
  #Return
	return(new.x)
}

##' @param object A \code{\link{comparative.comm}} object
##' @return Names of the traits or environmental variables
##' @rdname cc.manip
##' @export
trait.names <- function(object) {
    if(is.null(object$data)) return(NULL)
    colnames(object$data)
}

##' @export
##' @rdname cc.manip
env.names <- function(object) {
    if(is.null(object$env)) return(NULL)
    colnames(object$env)
}

##' @export
##' @rdname cc.manip
species <- function(x){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    return(colnames(x$comm))
}

##' @export
##' @rdname cc.manip
`species<-` <- function(x, value){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    colnames(x$comm) <- value
    x$phy$tip.label <- value
    if(!is.null(x$data))
        rownames(x$data) <- value
    return(x)
}
##' @export
##' @rdname cc.manip
sites <- function(x){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    return(rownames(x$comm))
}
##' @export
##' @rdname cc.manip
`sites<-` <- function(x, value){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    rownames(x$comm) <- value
    if(!is.null(x$env))
        rownames(x$env) <- value
    return(x)
}

##' @export
##' @rdname cc.manip
`traits<-` <- function(x, value){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    x$data <- value
    return(comparative.comm(x$phy, x$comm, x$data, x$env))
}

##' @export
##' @rdname cc.manip
traits <- function(x){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    return(x$data)
}

##' @export
##' @rdname cc.manip
`env<-` <- function(x, value){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    x$env <- value
    return(comparative.comm(x$phy, x$comm, x$data, x$env))
}

##' @export
##' @rdname cc.manip
env <- function(x){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    return(x$env)
}

##' @export
##' @rdname cc.manip
`comm<-` <- function(x, value){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    x$comm <- value
    return(comparative.comm(x$phy, x$comm, x$data, x$env))
}

##' @export
##' @rdname cc.manip
comm <- function(x){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    return(x$comm)
}

##' @export
##' @rdname cc.manip
tree <- function(x){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    return(x$phy)
}

##' @export
##' @rdname cc.manip
phy <- tree

##' @export
##' @rdname cc.manip
`tree<-` <- function(x, value){
    if(!inherits(x, "comparative.comm"))
        stop("'", deparse(substitute(x)), "' not of class 'comparative.comm'")
    x$phy <- value
    return(comparative.comm(x$phy, x$comm, x$data, x$env))
}

##' @export
##' @rdname cc.manip
`phy<-` <- `tree<-`

#' @param data A \code{\link{comparative.comm}} object
#' @export
#' @rdname cc.manip
assemblage.phylogenies <- function(data){
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    subtrees <- vector("list", nrow(data$comm))
    for(i in seq(nrow(data$comm)))
        subtrees[[i]] <- drop_tip(data$phy, rownames(data$comm)[data$comm[i,]!=0])
    return(subtrees)
}

#' @export
#' @rdname cc.manip
#' @param abundance.weighted whether to create to create a
#' @param row.names ignored
#' @param optional ignored
#' presence-absence dataset (default: FALSE)
as.data.frame.comparative.comm <- function(x, row.names=NULL, optional=FALSE, abundance.weighted=FALSE, ...){
    #Argument handling
    if(!inherits(x, "comparative.comm"))  stop("'data' must be a comparative community ecology object")

    #Wrapper for expanding
    # - tricky because of dummy factor variables
    expand <- function(data, env, n.spp, n.sites){
        if(!is.null(data)){
            mat <- sapply(data, function(x) model.matrix(~x-1))
            if(is.list(mat))
                mat <- do.call(cbind, mat)
            y <- 1
            for(i in seq(ncol(data))){
                if(is.character(data[,i]) | is.factor(data[,i])){
                    colnames(mat)[y:(y+length(unique(data[,i]))-1)] <- paste(names(data)[i], unique(data[,i]), sep=".")
                    y <- y+length(unique(data[,i]))
                } else {
                    colnames(mat)[y] <- names(data)[i]
                    y <- y+1
                }
            }
            if(env) mat <- apply(mat, 2, rep, each=n.spp) else mat <- apply(mat, 2, rep, n.sites)
        } else mat <- matrix(nrow=prod(n.spp, n.sites),ncol=0)
        return(mat)
    }
    
    #Make matrices - note that dummy variables make this tricky
    output <- matrix(t(x$comm), ncol=1)
    env <- expand(x$env, TRUE, length(species(x)), length(sites(x)))
    traits <- expand(x$data, FALSE, length(species(x)), length(sites(x)))
    site <- matrix(t(row(x$comm)), ncol=1)
    
    #Format output and return
    output <- as.data.frame(cbind(output, env, traits, site))
    rownames(output) <- NULL
    if(abundance.weighted){
        names(output)[1] <- "abundance"} else {
            names(output)[1] <- "presence"
            output[,1][output[,1] > 1] <- 1
        }
    names(output)[ncol(output)] <- "site"
    return(as.data.frame(output))
}
