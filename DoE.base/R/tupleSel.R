tupleSel <- function(design, type="complete", selprop=0.25, ...) UseMethod("tupleSel")
tupleSel.design <- function(design, type="complete", selprop=0.25, ...){
    ## eventually, there might be more possible types,
    ## automatic default type selection, ...
    if (!"design" %in% class(design)) stop("tupleSel.design is for class design objects only")
    tupleSel.default(design[,names(factor.names(design))], type=type, selprop=selprop, ...)
}
tupleSel.default <- function(design, type="complete", selprop=0.25, ...){
    ## so far, types complete, worst.abs and worst.rel implemented
    ##
    hilf <- round(GWLP(design), 10)
    dim <- min(which(hilf[-1] > 0))
    if (dim <= 2) stop("the design is not orthogonal, \noption", type, "is not available")
    if (dim > 4) stop("the design has resolution larger than IV, \noption", type, "is not available")

    if (dim==3) fn <- "P3.3" else fn <- "P4.4"
    rela <- (type %in% c("complete", "worst.rel"))
    parft <- type=="worst.parft"
    parftdf <- type=="worst.parftdf"
    hilf.all <- eval(parse(text=paste(fn, "(design, rela=",rela,", parft=",parft,", parftdf=",parftdf,", detailed=TRUE)")))
    hilf <- attr(hilf.all, "detail")
    attr(hilf.all, "detail") <- NULL
    if (type=="complete") pos <- which(hilf==1)
    else {
       pos <- which(hilf > quantile(hilf, 1-selprop))
       pos <- pos[ord(data.frame(hilf[pos]), decreasing=TRUE)]
       }
    aus <- lapply(strsplit(names(hilf[pos]), ":"), as.numeric)
    if (rela) attr(aus, "RPFT") <- hilf.all else 
    if (parft) attr(aus, "PARFT") <- hilf.all else 
    if (parftdf) attr(aus, "PARFTdf") <- hilf.all else 
    attr(aus, "PFT") <- hilf.all
    attr(aus, "detail") <- hilf[pos]
    aus
}

