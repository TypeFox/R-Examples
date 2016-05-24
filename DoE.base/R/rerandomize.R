rerandomize.design <- function(design, seed=NULL, block=NULL, ...){
    ## function to re-randomize a design object
    ## particularly interesting for replicated designs,
    ##    if users don't want to have the randomization in blocks
    ## or for blocked designs, if users want the blocks in randomized order
    ## or for designs for which a block factor is included as an experimental 
    ##    factor and has not been randomized upon
    if (!"design" %in% class(design))
        stop("the function works on class design objects only")
    di <- design.info(design)
    if (!is.null(di$response))
        stop("the design has responses already and must not be re-randomized.")
    if (!is.null(seed)){
       if (!is.numeric(seed)) stop("seed must be a number")
       if (!length(seed)==1) stop("seed must be a single number")
       }
    if (!is.null(block)){
       if (!is.character(block)) 
       stop("block must be a character string (name of the block factor)")
       if (length(block)>1) stop("block must yield a single factor name")
       bf <- which(names(di$factor.names)==block)
       if (length(bf)==0) stop("block does not specify an experimental factor")
    }
    ro <- run.order(design)
    desnum <- desnum(design)
    design <- undesign(design)

    aw <- FALSE
    awro <- FALSE
    bl <- FALSE
    sp <- FALSE
    if (di$replications>1 & !di$repeat.only) aw <- TRUE
    if (di$replications>1 & di$repeat.only) awro <- TRUE
    if (length(grep("blocked", di$type, fixed=TRUE))>0) bl <- TRUE
    if (length(grep("splitplot", di$type, fixed=TRUE))>0) sp <- TRUE
    if ((bl || sp || aw || awro) && !is.null(block)) 
       stop("block can only be specified for unblocked and unreplicated designs without split plot structure")

    if (!is.null(seed)) set.seed(seed)
    ## replications without complications
    if (!(bl | sp | !is.null(block))){
      rp <- di$replications
      if (awro){
        nr <- nrow(design)%/%rp
        neworder <- sample(nr)
        neworder <- rep((neworder-1)*rp, each=rp) + rep(1:rp, nr)
        }
      else{
        nr <- nrow(design)
        neworder <- sample(nr)
        }
      }
    else {
      ## splitplot
      if (sp) {
        ps <- di$plotsize
        rp <- di$replications
        if (!awro){
          # plot <- as.numeric(getblock(design, combine=TRUE))
          nr <- nrow(design)%/%ps
          neworder <- sample(nr)
          neworder <- rep((neworder-1)*ps, each=ps) + unlist(lapply(1:nr, function(obj) sample(ps)))
        }
        else{
          # plot <- as.numeric(getblock(design)$plots)
          nr <- nrow(design)%/%(ps*di$replications)
          neworder <- sample(nr)
          neworder <- rep((neworder-1)*ps, each=ps*rp) +
               rep(unlist(lapply(1:nr, function(obj) sample(ps))), each=rp) +
               rep(1:rp, nr*ps)
        }
      }  ## end of split plot
    if (bl) {
         ## block
        bs <- di$blocksize
        nb <- di$bbreps
        nw <- di$wbreps
        if (!awro) {
          nr <- nrow(design)%/%(bs*nw)
          neworder <- sample(nr)
          neworder <- rep((neworder-1)*bs*nw, each=bs*nw) + unlist(lapply(1:nr, function(obj) sample(bs*nw)))
        }
        else{
          nr <- nrow(design)%/%(bs*nw)
          neworder <- sample(nr)
          neworder <- rep((neworder-1)*bs, each=bs*nw) +
               rep(unlist(lapply(1:nr, function(obj) sample(bs))), each=nw) +
               rep(1:nw, nr*bs)
        }
    }  ## end of block
    if (!is.null(block)){
        di$block.name <- block
        di$nblocks <- nr <- length(di$factor.names[[bf]])
        di$blocksize <- bs <- di$nruns%/%di$nblocks
        di$nfactors <- di$nfactors-1
        di$factor.names <- di$factor.names[-bf]
        di$bbreps <- 1
        di$wbreps <- 1
        if (!is.null(di$nlevels)) di$nlevels <- di$nlevels[-bf]
        if (!is.null(di$selected.columns)) di$selected.columns <- di$selected.columns[-bf]
        di$type <- paste(di$type, "blocked", sep=".")
        stdorder <- ord(ro)  ## brings design to standard order (so far unblocked should be unique)
        blockorder <- ord(design[stdorder,bf, drop=FALSE]) ## ordered by block, not yet randomized
        neworder <- sample(nr)
        neworder <- rep((neworder-1)*bs, each=bs) + unlist(lapply(1:nr, function(obj) sample(bs)))
        ## combine stdorder with blockorder and neworder
        neworder <- stdorder[blockorder][neworder]
        design <- cbind(design[,bf],design[,-bf])
        names(design)[1] <- block
        contrasts(design[[block]]) <- "contr.treatment"
        di$creator <- list(di$creator, paste("rerandomized with block factor", block, "and seed", seed))
    }  ## end of experimental block factor                                             
    }  ## end of block or splitplot or experimental block factor
   design <- design[neworder,]
   ro <- ro[neworder,]
   desnum <- desnum[neworder,]
   if (!is.null(block)) {
      hilf <- as.numeric(as.character(ro[[1]]))
      hilfrank <- unlist(lapply(1:nr, function(obj) rank(hilf[(obj-1)*bs+(1:bs)])))
      ro[[1]] <- ro[[3]] <- factor(paste(hilf,as.numeric(as.factor(design[[1]])),hilfrank,sep="."))
      desnum <- model.matrix(~., design)
      ro$run.no <- 1:nrow(ro)
      row.names(design) <- ro$run.no
      row.names(ro) <- ro$run.no
      }
   attr(design, "desnum") <- desnum
   attr(design, "run.order") <- ro
   attr(design, "design.info") <- di
   class(design) <- c("design", "data.frame")
   design
}