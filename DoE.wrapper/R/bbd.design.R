bbd.design <- function(nfactors, ncenter=4, factor.names = NULL, default.levels=c(-1,1), 
          block.name=NULL, randomize=TRUE, seed=NULL, ...){
   creator <- sys.call()
    ##  use factor.names for creating the coding option coding=make.formulas(paste("x",1:nfactors,sep=""),factor.names)
    ## incorporate checks on factor.names from lhs and randomization / replication from FrF2 or pb
    ## output the design in the format of all other packages
    ## bbd(nfactors, n0=ncenter, randomize=randomize, coding=make.formulas(paste("x",1:nfactors,sep=""),factor.names))
    if (!nfactors %in% 3:7) stop("Box-Behnken designs are implemented for 3 to 7 factors only.")
    if (!is.numeric(ncenter)) stop("ncenter must be a number.")
    if (!floor(ncenter)==ncenter) stop("ncenter must be an integer number.")
    if (!is.numeric(default.levels)) stop("default.levels must be numeric")
    if (!length(default.levels)==2) stop("there must be 2 default.levels (the center is calculated)")
    if (!length(unique(default.levels))==2) stop("duplicate default levels")
    if (!(is.character(block.name) | is.null(block.name))) stop("block.name must be NULL or character")
    if (!is.logical(randomize)) stop("randomize must be logical")
    if (!(nfactors %in% c(4,5) | is.null(block.name))) {
         warning("Box-Behnken designs can only be blocked in case of 4 or 5 factors. No blocking was done.")
         block.name <- NULL
         }
         if (!is.null(block.name)) block.name <- make.names(block.name)
    if (is.character(factor.names)){
        if (!length(factor.names)==nfactors) stop("mismatch between nfactors and factor.names")
        if (!all(unique(factor.names)==factor.names)) stop("duplicate factor names")
        hilf <- rep(list(default.levels),nfactors)
        names(hilf) <- factor.names
        factor.names <- hilf
        }
    if (is.null(factor.names)){ 
      factor.names <- rep(list(default.levels),nfactors)
      names(factor.names) <- Letters[1:nfactors]
    }
    if (!is.list(factor.names)) stop("if given, factor.names must be a character vector or a list")
    if (!length(factor.names)==nfactors) stop("mismatch between nfactors and factor.names")
    for (i in 1:nfactors) if (identical(factor.names[[i]],"")) factor.names[[i]] <- default.levels
    if (is.list(factor.names) & !length(unique(names(factor.names)))==nfactors)
            names(factor.names) <- Letters[1:nfactors]
    ## make all factor names valid R names
    names(factor.names) <- make.names(names(factor.names), unique=TRUE)

    if (randomize & !is.null(seed)) set.seed(seed)
    aus <- .bbd.1.41(nfactors, n0=ncenter, 
        block = if (is.null(block.name)) FALSE else block.name, randomize=randomize, 
        coding=make.formulas(paste("x",1:nfactors,sep=""),factor.names))
    ## must still be made into design with the usual information available
    design <- decode.data(aus)
    class(design) <- c("design","data.frame")
    desnum <- aus
    colnames(desnum) <- colnames(design)
    desnum <- model.matrix(~.,desnum)[,-1]
    rownames(design) <- rownames(desnum) <- 1:nrow(aus)
    desnum(design) <- desnum
    ## changes added 27 01 2011 (for proper ordering of design)
    orig.no.levord <- sort(as.numeric(rownames(aus)),index=TRUE)$ix
    orig.no <- factor(rownames(aus), levels=unique(rownames(aus)[orig.no.levord]))
    run.order(design) <- data.frame(run.no.in.std.order=orig.no, run.no=1:nrow(aus), run.no.std.rp=rownames(aus))
    di <- list(type="bbd", nruns=nrow(design), nfactors=nfactors, factor.names=factor.names, quantitative=rep(TRUE, nfactors))
    if (!is.null(block.name)){ blocklist <- list(nblocks=if(nfactors==4) 3 else 2,
               blocksize=if(nfactors==4) 12 else 24, block.name=block.name, bbreps=1, wbreps=1)
          di$type <- "bbd.blocked"
          di <- c(di,blocklist)
          }
    di <- c(di, list(randomize=randomize, seed=seed, replications=1, repeat.only=FALSE, ncenter=ncenter, creator=creator, 
        coding=make.formulas(paste("x",1:nfactors,sep=""),factor.names)))
    names(di$coding) <- paste("x",1:nfactors,sep="")
    design.info(design) <- di
    design
}
