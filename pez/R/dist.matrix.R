#' Make co-existence matrices based on phylogeny (and/or) traits, and
#' community or environemntal overlap
#' 
#' @param x an object
#' @param ... not used
#' @details \code{comm.dist} returns the 1 - co-existence of
#' species. Look at how this is calcualted; it incorporates
#' abundances, and if you don't want it to do so simply call it on a
#' presence/absensence (1/0) matrix.
#' @rdname dist.xxx
#' @name dist.xxx
#' @export
comm.dist <- function(x) UseMethod("comm.dist", x)
#' @export
#' @rdname dist.xxx
#' @importFrom stats as.dist
comm.dist.matrix <- function(x){
	output <- matrix(ncol=ncol(x), nrow=ncol(x))
	for(i in seq(ncol(x))){
		for(j in seq(ncol(x))){
			output[i,j] <- 1 - sum(x[,i]>0 & x[,j]>0) / sum(x[,i]>0 | x[,j]>0)
		}
	}
        if(any(!is.finite(output))){
            output[!is.finite(output)] <- 1
            warning("Co-existence matrix contains non-overlapping species; ",
                    "treating as maximally dissimilar")
        }
        return(as.dist(output))
}
#' @export
#' @rdname dist.xxx
comm.dist.comparative.comm <- function(x) return(comm.dist(x$comm))

#' @details \code{traits.dist} returns the functional trait distance
#' of species
#' @param dist.func a function for computing distances.  The default,
#' \code{dist.func.default}, returns a Euclidean distance of the
#' scaled and centred data.
#' @param alltogether should one multivariate distance matrix be
#' computed for all traits at once (DEFAULT; \code{alltogether =
#' TRUE}) or for each trait at a time (\code{alltogether = FALSE})?
#' @rdname dist.xxx
#' @export
traits.dist <- function(x, dist.func = dist.func.default, ...) UseMethod("traits.dist", x)
#' @rdname dist.xxx
#' @export
traits.dist.comparative.comm <- function(x, dist.func = dist.func.default, alltogether = TRUE, ...){
    if(is.null(x$data)) stop("No trait data for which to compute a trait distance matrix")
    if(alltogether){
        return(traits.dist(x$data))
    } else {
        return(sapply(as.data.frame(x$data), dist.func)) 
    }
}
#' @export
#' @rdname dist.xxx
traits.dist.default <- function(x, dist.func = dist.func.default, ...) dist.func(x)
#' @export
#' @rdname dist.xxx
traits.dist.data.frame <- function(x, dist.func = dist.func.default, ...) dist.func(as.matrix(x))
#' @export
#' @rdname dist.xxx
#' @importFrom stats dist
dist.func.default <- function(x) dist(scale(x, center=TRUE, scale=TRUE))

#' @details \code{phylo.dist} returns the phylogenetic (cophenetic)
#' distance of species
#' @export
#' @rdname dist.xxx
phylo.dist <- function(x, ...) UseMethod("phylo.dist", x)
#' @export
#' @rdname dist.xxx
#' @importFrom stats as.dist
#' @importFrom ape cophenetic.phylo
phylo.dist.phylo <- function(x, ...) as.dist(cophenetic.phylo(x))
#' @export
#' @rdname dist.xxx
phylo.dist.comparative.comm <- function(x, ...) phylo.dist(x$phy)


#' @details \code{funct.phylo.dist} returns the combined phylogenetic
#' and trait distances of species, based on the traitgram approach of
#' Cadotte et al. (2013).
#' @param phyloWeight phylogenetic weighting parameter (referred to as
#' \code{a} in Cadotte et al. (2013)
#' @param p exponent giving the exponent to use for combining
#' functional and phylogenetic distances (the default, \code{p = 2},
#' gives a Euclidean combination).
#' @details Make functional phylogenetic distance matrix
#' @references Cadotte M.A., Albert C.H., & Walker S.C. The ecology of differences: assessing community assembly with trait and evolutionary distances. Ecology Letters 16(10): 1234--1244.
#' @rdname dist.xxx
#' @export
funct.phylo.dist <- function(x, phyloWeight, p = 2, ...) UseMethod("funct.phylo.dist", x)
#' @export
funct.phylo.dist.comparative.comm <- function(x, phyloWeight, p, ...) {
    #Assertions and argument handling
    if(!inherits(x, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    if(is.null(x$data)) stop("'data' must contain trait data")
    if(phyloWeight < 0 | phyloWeight > 1)
        stop("'phyloWeight' must be between 0 and 1")
    if(!is.numeric(p)) stop("'p' must be a numeric")
    
    FDist <- traits.dist(x)
    PDist <- phylo.dist(x)
    FDist <- FDist/max(FDist)
    PDist <- PDist/max(PDist)
    (phyloWeight * PDist^p + (1 - phyloWeight) * FDist^p)^(1/p)
}

#' @details \code{pianka.dist} returns the environemntal tolerances
#' distance matrices of species. Based on Pianka's distance (i.e.,
#' niche overlap based on environmental variables at co-occuring
#' sites), as defined in Cavender-Bares et al. (2004) - likely not the
#' original reference!
#' @export
#' @rdname dist.xxx
#' @references Cavender-Bares J., Ackerly D.D., Baum D.A. & Bazzaz F.A. (2004) Phylogenetic overdispersion in Floridian oak communities. The Americant Naturalist 163(6): 823--843.
pianka.dist <- function(x, ...) UseMethod("pianka.dist", x)

#' @rdname dist.xxx
#' @param env environmental variable to be used to calculate the
#' distance matrix
#' @export
#' @importFrom stats as.dist
pianka.dist.matrix <- function(x, env = NULL, ...){
    #Checks and assertions
    if(!is.numeric(x)) stop("Need numeric community matrix for Pianaka calculations")
    if(!is.factor(env)) stop("Pianaka calculations require a factor as the second argument")
    #Find the proportional occupancy
    propOcc <- matrix(nrow=ncol(x), ncol=length(levels(env)))
    colnames(propOcc) <- levels(env)
    #Matrices with one row become vectors; this confuses colSums
    for(j in seq(ncol(propOcc)))
        if(sum(env==colnames(propOcc)[j])>1) propOcc[,j] <- colSums(x[env==colnames(propOcc)[j],]) else propOcc[,j] <- x[env==colnames(propOcc)[j],]
    propOcc <- apply(propOcc, 2, function(x) x/rowSums(propOcc))
    #Get the Pianka overlap
    pianka <- matrix(ncol=ncol(x), nrow=ncol(x))
    #Matrix is symmetrical, so use this to speed things up
    for(j in seq(from=1, to=ncol(x)-1)){
        for(k in seq(from=j+1, to=ncol(x))){
            pianka[j,k] <- sum(propOcc[j,] * propOcc[k,]) / sqrt(sum(propOcc[j,]^2) + sum(propOcc[k,]^2))
            pianka[k,j] <- pianka[j,k]
        }
    }
    return(as.dist(pianka))
}

#' @export
#' @rdname dist.xxx
#' @importFrom stats as.dist
pianka.dist.comparative.comm <- function(x, alltogether = TRUE, ...){
	#Checks and assertions
	if(is.null(x$env)) stop("Cannot calculate Pianka distances without environmental data")
	if(any(!sapply(x$env, is.factor))) stop("Cannot calculate Pianka distances of environmental data non-factor-level environmental data")
	#Pianaka matrix for each environmental variable
	pianka <- array(dim=c(ncol(x$comm), ncol(x$comm), ncol(x$env)))
	for(i in seq(ncol(x$env)))
		pianka[,,i] <- as.matrix(pianka.dist.matrix(x$comm, x$env[,i]))
  if(alltogether)
    return(as.dist(apply(pianka, 1:2, mean))) else return(pianka)
}
