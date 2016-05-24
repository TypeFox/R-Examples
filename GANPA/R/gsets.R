preproc.gsets <-
function (gsets, gbg = NULL, size.max = 500, size.min = 15) 
{
    if (!is.null(gbg)) {
        gbg <- unique(gbg)
        gsets <- lapply(gsets, function(x) {
            intersect(x, gbg)
        })
    }
    else {
        gsets <- lapply(gsets, unique)
    }
    size <- unlist(lapply(gsets, length))
    gsets[which(size <= size.max & size >= size.min)]
}
weight.gset.options <-
function (isets, gset, glist = NULL, option) 
{
    gcatalog <- names(isets)
    N <- length(gcatalog)
    if (is.null(glist)) 
        glist <- gset
    gset <- gset[!is.na(match(gset, gcatalog))]
    M <- length(gset)
    glist.ind <- match(glist, gcatalog)
    w <- rep(NA, length(glist.ind))
    for (i in 1:length(glist.ind)) {
        if (is.na(glist.ind[i])) 
            next
        gi <- isets[[glist.ind[i]]]
        K <- length(gi)
        X <- length(intersect(gi, gset))
        X0 <- M * K/N
        w[i] <- switch(option, x = X, x0 = X0, dx = X - X0, or = X/X0, 
            stop("Unknown 'option'"))
    }
    names(w) <- glist
    return(w)
}
weight.gset.test <-
function (isets, gset, glist = NULL) 
{
    dx <- weight.gset.options(isets = isets, gset = gset, glist = glist, 
        option = "dx")
    dx[which(dx < 0 | is.na(dx))] <- 0
    log2(dx + 2)
}
weight.gsets.test <-
function (isets, gsets) 
{
    w <- lapply(gsets, function(gset) {
        weight.gset.test(isets = isets, gset = gset)
    })
}
read.gmt <-
function (file.gmt) 
{
    gmt <- read.table(file.gmt, header = FALSE, sep = "\n", as.is = TRUE, 
        quote = "\"", comment.char = "")[, 1]
    gmt <- strsplit(gmt, split = "\t")
    id <- unlist(lapply(gmt, function(gs) {
        gs[1]
    }))
    name <- unlist(lapply(gmt, function(gs) {
        gs[2]
    }))
    GSEA2symbols <- lapply(gmt, function(gs) {
        sort(unique(gs[-(1:2)]))
    })
    names(GSEA2symbols) <- id
    GSEA_id2name <- name
    names(GSEA_id2name) <- id
    return(list(id2gene = GSEA2symbols, id2name = GSEA_id2name))
}
write.gmt <-
function (gs.list, file) 
{
    if (length(grep(".gmt$", file)) == 0) 
        stop("write.gmt(): file should be .gmt")
    gs.names <- names(gs.list)
    gs.desc <- rep(NA, length(gs.names))
    gs.lines <- unlist(lapply(gs.list, function(x) {
        paste(x, collapse = "\t")
    }))
    gs.lines <- paste(gs.names, gs.desc, gs.lines, sep = "\t")
    writeLines(gs.lines, con = file)
}
write.gwt <-
function (gs.weight, file) 
{
    if (length(grep(".gwt$", file)) == 0) 
        stop("write.gwt(): file should be .gwt")
    gs.names <- names(gs.weight)
    gs.desc <- rep(NA, length(gs.names))
    gs.lines <- unlist(lapply(gs.weight, function(x) {
        paste(names(x), round(x, 3), sep = "|", collapse = "\t")
    }))
    gs.lines <- paste(gs.names, gs.desc, gs.lines, sep = "\t")
    writeLines(gs.lines, con = file)
}
