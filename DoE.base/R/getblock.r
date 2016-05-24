getblock <- function(design, combine=FALSE, ...){
    di <- design.info(design)
    aw <- FALSE
    awro <- FALSE
    bl <- FALSE
    sp <- FALSE
    if (di$replications>1 & !di$repeat.only) aw <- TRUE
    if (di$replications>1 & di$repeat.only) awro <- TRUE
    if (length(grep("blocked", di$type, fixed=TRUE))>0) bl <- TRUE
    if (length(grep("splitplot", di$type, fixed=TRUE))>0) sp <- TRUE
    if (bl) {
            if (di$bbreps>1 | (di$wbreps>1 & !di$repeat.only)) aw <- TRUE
            if (di$wbreps>1 & di$repeat.only) awro <- TRUE
            }
    if (bl & !(aw | awro)) stop("Nothing was done, the design contains the appropriate block factor ", di$block.name)
    if ((!aw) & (!awro) & (!sp)) stop("Nothing was done, as the design does not contain replications or repeated measurements")
    if (combine & !(bl | (sp & (aw | awro)))) combine <- FALSE  ## nothing to combine for aw only
    ro <- run.order(design)
    if (!all(ro$run.no == sort(ro$run.no))) 
        stop("getblock does not work for designs that have been reordered after creation")
    rov <- ro$run.no.std.rp
    if (is.factor(rov)) rov <- as.character(rov)
    rovs <- strsplit(rov,".",fixed=TRUE)
    rovs <- lapply(rovs, as.numeric)
    
    ncenter <- di$ncenter
    if (is.null(di$ncenter)) ncenter <- 0

    ## replications only
    if (!(bl | sp)) {
        if (aw) blocks <- as.factor(sapply(rovs, function(obj) obj[2]))
        if (awro) blocks <- as.factor(rep(1:(di$nruns+ncenter), each=di$replications)) 
        }
    else{
       if (bl){
         ## blocked designs
         ## distinguish within and between or both
         blocks <- factor(sapply(rovs, function(obj) obj[2]), levels=1:di$nblocks)
         if (di$bbreps > 1) between.reps <- factor(sapply(rovs, function(obj) obj[4]), levels=1:di$bbreps)
         else { ## then di$wbreps must be larger than 1
            if (aw) within.reps <- factor(sapply(rovs, function(obj) obj[4]), levels=1:di$wbreps)
            if (awro) within.reps <- as.factor(rep(1:(di$blocksize+ncenter), each=di$wbreps,times=di$nblocks*di$bbreps)) 
         }
         ## both
         if (di$bbreps > 1 & di$wbreps > 1){ 
            if (aw) within.reps <- factor(sapply(rovs, function(obj) obj[5]), levels=1:di$wbreps)
            if (awro) within.reps <- as.factor(rep(1:(di$blocksize+ncenter), each=di$wbreps,times=di$nblocks*di$bbreps)) 
            }
         blocks <- data.frame(blocks)
         if (exists("between.reps")) blocks <- cbind(blocks, between.reps)
         if (exists("within.reps")) blocks <- cbind(blocks, within.reps)
       }
       if (sp){
         ## split plot
         blocks <- factor(sapply(rovs, function(obj) obj[2]), levels=1:di$nWPs)
         if (aw)
         blocks <- data.frame(plots=blocks, reps=sapply(strsplit(rov,".",fixed=TRUE), function(obj) obj[4]))
         if (awro)
         blocks <- data.frame(plots=blocks, reps=as.factor(rep(1:di$plotsize, each=di$replications,times=di$nWPs)))
      }              
      if (combine) {
            ## create single factor from multi-column data frame
            reihenfolge <- ord(blocks)
            blocks <- apply(as.matrix(blocks), 1, function(obj) paste(obj,collapse="."))
            blocks <- factor(blocks, levels=unique(blocks[reihenfolge]))
          }
    }
    blocks
}
