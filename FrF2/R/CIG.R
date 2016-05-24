CIG <- function (design, select.catlg = catlg, static = FALSE, layout = layout.auto,
    label = "num", plot = TRUE, ...)
{

    ## updates needed for documentation:
    ##         add requirement set examples

    ## permit handover of vertex labels
    ll <- list(...)
    if ("vertex.label" %in% names(ll)) vertex.label <- ll$vertex.label

    ## create graph picture for design entry
    if ("catlg" %in% class(design)) {
        if (length(design) > 1)
            stop("design must not contain more than one catlg entry")
        design <- design[[1]]
    }
    else {
        if (is.matrix(design) || (is.character(design) & length(design)>1) || is(design,"formula")){
        ## Formel
        if (is(design,"formula")){
            fn <- row.names(attr(terms(formula(design)), "factors"))
            design <- estimable.check(design, length(fn), fn)
            names(design) <- c("clear.2fis","nfac")
             if (!exists("vertex.label", inherits=FALSE))  vertex.label <- fn
             else if (!length(vertex.label)==length(fn)) warning("vertex.label has wrong length")
            design$res <- 4
        } 
        ## character or numeric two-row matrix, or character vector with elements from Letters
        if (is.matrix(design)){
             if (!nrow(design)==2) stop("matrix design must have two rows.")
             if (any (design[1,]==design[2,])) stop("entries in the same column of matrix design must be different")
             fn <- names(table(design))
             if (!exists("vertex.label", inherits=FALSE))  vertex.label <- fn
             else if (!length(vertex.label)==length(fn)) warning("vertex.label has wrong length")
             design <- list(clear.2fis=matrix(sapply(design, function(obj) which(fn==obj)),nrow=2), nfac=length(fn),res=4)
             }
        if (is.character(design) && length(design)>1){
             if (!all(nchar(design)==2)) stop("character vector design must have length 2 entries only")
             design <- estimable.check(design,NULL,NULL)
             if (!exists("vertex.label", inherits=FALSE))  vertex.label <- Letters[1:design$nfac]
             else if (!length(vertex.label)==design$nfac) warning("vertex.label has wrong length")
             names(design) <- c("clear.2fis","nfac")
             if (any (design$clear.2fis[1,]==design$clear.2fis[2,]))
               stop("characters in the same element of vector design must be different")
             design$res <- 4
             }
             ## now, design is a two-row numeric matrix with entries from 1 to nfac
        }
    else{
        if (!"catlg" %in% class(select.catlg))
            stop("select.catlg must be a catalogue")
        if (!(is.character(design) & length(design) == 1))
            stop("design must be a design name")
        if (!(design %in% names(select.catlg)))
            stop("design must be a design name that occurs in select.catlg")
        design <- select.catlg[[design]]
    }
    }
    if (!exists("vertex.label", inherits=FALSE)) {
        vertex.label <- 1:design$nfac
        if (!label == "num")
            vertex.label <- Letters[vertex.label]
    }
    go2 <- graph.empty(n = design$nfac, directed = FALSE)
    if (!length(design$clear.2fis) == 0)
        go2 <- add.edges(go2, design$clear.2fis)
    if (design$res < 4)
        warning("the design is of resolution less than IV")
    if (plot) {
        if (!static) {
            id <- tkplot(go2, vertex.label = vertex.label, ...)
            invisible(list(graph = go2, coords = tkplot.getcoords(id)))
        }
        else {
            invisible(go2)
            plot(go2, layout = layout, vertex.label = vertex.label,
                ...)
        }
    }
    else return(go2)
}

CIGstatic <- function(graph, id, label = "num", xlim = c(-1, 1), ylim = c(1, 
    -1), ...){
    ## get coordinates for static graph from dynamic one
    if ("list" %in% class(graph)) 
        if (names(graph)[1] == "graph") graph <- graph$graph
    if (!exists("vertex.label", inherits = FALSE)) {
        vertex.label <- 1:graph[[1]][[1]][2]
        ## was 1:graph[[1]], changed 18 Apr 2013
        if (!label == "num") 
            vertex.label <- Letters[vertex.label]
    }
    coords <- tkplot.getcoords(id)
    plot(graph, layout = coords, vertex.label = vertex.label, 
        xlim = xlim, ylim = ylim, ...)
}
