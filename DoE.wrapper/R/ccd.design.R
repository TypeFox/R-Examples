ccd.design <- function(nfactors=NULL, factor.names=NULL, default.levels=c(-1,1), ncube=NULL, 
        resolution=if (identical(blocks,1) & is.null(ncube)) 5 else NULL, 
        generators=NULL, ncenter = 4, alpha = "orthogonal", 
        replications=1, block.name="Block.ccd", blocks=1, 
        randomize=TRUE, seed=NULL, ...){
    ## can make n.c omittable
    ## by defaulting to 2^nfactors or better to the smallest number that permits a resolution V design
    
    ## use FrF2 for defining the generators (ideally select a resolution V design)
    ## use log2(n.c) for basis
    ##  use factor.names for creating the coding option coding=make.formulas(paste("x",1:nfactors,sep=""),factor.names)

    creator <- sys.call()
    if (!is.null(ncube)) if (!is.numeric(ncube)) stop("ncube must be numeric")
    if (!is.null(ncube)){ 
         k <- round(log2(ncube))
         if (!2^k==ncube) stop("ncube must be a power of 2")
        }
    if (is.null(ncube) & blocks>1) stop("blocks>1 requires specification of ncube")
    if (is.null(nfactors) & is.null(factor.names)) stop("nfactors or factor.names must be given")
    if (is.null(nfactors)){
         if (!(is.character(factor.names) | is.list(factor.names)))
             stop("factor.names must be a character vector or a list")
             nfactors <- length(factor.names)
         }
    else if (is.null(factor.names)) {
                 factor.names <- rep(list(default.levels),nfactors)
                 names(factor.names) <- paste("X",1:nfactors,sep="")
         }
    block.name <- make.names(block.name) ## make block name a valid R name
    #if (blocks>1 & bbreps*wbreps>1) stop("replicated blocked designs are currently not covered by ccd.design")
    if (is.list(factor.names) & !length(unique(names(factor.names)))==nfactors)
            names(factor.names) <- paste("X",1:nfactors,sep="")
    if (!length(factor.names)==nfactors) stop("mismatch between nfactors and length of factor.names")
    ## make all factor names valid R names
    names(factor.names) <- make.names(names(factor.names),unique=TRUE)

    if (!is.numeric(ncenter)) stop("ncenter must be numeric")
    if (!all(ncenter==floor(ncenter))) stop("ncenter must be integer")
    if (!length(ncenter) %in% c(1,2)) stop("ncenter must have one or two elements")
    
    if (!(is.null(seed) | is.numeric(seed))) stop("seed must be NULL or numeric")
    if (!is.null(seed)){
      if (!all(floor(seed)==seed)) stop("non-NULL seeds must be integer")
      if (!length(seed)==2) stop("ccd.design needs a seed of length 2")
    }
    
    if (!(is.null(resolution) | is.null(generators))) 
        stop("resolution and generators cannot be specified together")
    if (!(is.null(resolution) | identical(blocks, 1))) 
        stop("resolution and blocks cannot be specified together")
    cube <- FrF2(nruns=ncube, nfactors=nfactors, factor.names=factor.names, 
          default.levels=default.levels, resolution=resolution, generators=generators, 
          block.name=block.name, blocks=blocks, 
          replications=replications, randomize = randomize, seed=seed[1])
    aus <- ccd.augment(cube, ncenter=ncenter, block.name=block.name, alpha=alpha, 
        randomize=randomize, seed=seed[2])
    hilf <- design.info(aus)
    hilf$creator <- creator
    design.info(aus) <- hilf
    aus
}

