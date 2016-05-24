iscube <- function (design, ...)
{
    ## TRUE, if cube points
    ## works for pb, FrF2 and ccd
    ## cube points do not include center points of cube portion (ccd)!
    if (!"design" %in% class(design))
        stop("this function is applicable for class design objects only")
    di <- design.info(design)
    if (!(length(grep("center", di$type)) > 0 | length(grep("ccd", di$type)) > 0))
        stop("this function requires a design with center points or a central composite design")
    ## determine center point positions
    ## modification Jan 2012 for making new function work on old design versions
    hilf <- as.character(run.order(design)$run.no.in.std.order)
    hilf2 <- rownames(design)
    ## handle old and new ccd designs (old ccd designs: row names are correct, 
    ##        run.order contains different row names)
    ## old ccd designs can only be handled in case of no replications
## new ccd designs (since 1.6-5): run.order is correct
    if (isTRUE(all.equal(sort(hilf),sort(hilf2)))) 
      if (!isTRUE(all.equal(hilf, hilf2))) hilf <- hilf2
    ## versions up ro FrF2 1.6-4 had run.no.in.std.order=0 for all center points
    ##    for later versions, the 0 may be accompanied by block or replication info
    aus <- (hilf == "0" | sapply(hilf, function(obj) unlist(strsplit(obj,".",fixed=TRUE))[1])=="0") | 
        (sapply(hilf, function (obj) length(grep("C", obj)) > 0) & 
          as.numeric(sapply(hilf, function(obj) substr(obj,4,99))) > di$ncube) |
          sapply(hilf, function (obj) length(grep("S", obj)) > 0)
    if (!sum(!aus) == di$ncube * di$replications[1]) {
        wrong <- TRUE
        if (!is.null(di$blocksize))
            if (sum(!aus) == di$ncube * di$bbreps[1] * di$wbreps[1] *
                di$nblocks[1])
                wrong <- FALSE
        if (wrong)
            stop("There is something wrong with the number of center points for this design.")
    }
    ## output cube point positions
    !aus
}

isstar <- function (design, ...)
{
    if (!"design" %in% class(design))
        stop("this function is applicable for class design objects only")
    di <- design.info(design)
    if (!length(grep("ccd", di$type)) > 0)
        stop("this function requires a central composite design")
    ## determine star point positions
    hilf <- as.character(run.order(design)$run.no.in.std.order)
    sapply(hilf, function (obj) length(grep("S", obj)) > 0) 
}

pickcube <- function(design, ...)
{
    aufruf <- sys.call()
    if (!"design" %in% class(design))
        stop("this function is applicable for class design objects only")
    di <- design.info(design)
    if (!length(grep("ccd", di$type)) > 0)
        stop("this function requires a central composite design")
    ro <- run.order(design)[!isstar(design),]
    desnum <- desnum(design)[!isstar(design),]
    hilf <- design[!isstar(design),]
    if (di$ncenter[1]>0) di$type <- "FrF2.center" else di$type <- "FrF2"
    di$nruns <- nrow(hilf)
    di$ncenter <- di$ncenter[1]
    di$nstar <- NULL
    di$creator <- append(di$creator, list(aufruf))
    attr(hilf,"design.info") <- di
    attr(hilf,"run.order") <- ro
    attr(hilf,"desnum") <- desnum
    class(hilf) <- c("design","data.frame")
    hilf
}