image.vhica <-
function (x, element = "", H1.test = "bilat", treefile = NULL, 
    skip.void = FALSE, species = NULL, p.threshold = 0.05, p.adjust.method = "bonferroni", 
    ncolors = 1024, main = element, threshcol = 0.1, colsqueeze=1, ...) 
{
    op <- par(no.readonly = TRUE)
    tree <- .prepare.phylo(treefile)
    if (!"phylo" %in% class(tree)) {
        species <- .check.species(x, user.species = species)
        tree <- list(tip.label = species)
    }
    else {
        species <- .check.species(x, user.species = species, 
            tree.species = tree$tip.label)
    }
    elements <- .element.present(x, element, species = tree$tip.label, 
        skip.void = skip.void)
    if(is.null(elements)) {
		stop("Element ", element, " cannot be found. Nothing to plot.")
    }
    if (length(tree) > 1 && skip.void) {
        missing.species <- tree$tip.label[!tree$tip.label %in% 
            sapply(strsplit(elements, ".", fixed = TRUE), function(el) el[1])]
        tree <- ape::drop.tip(tree, missing.species)
        species <- species[species %in% tree$tip.label]
    }
    stats <- .stat.matrix(vhica.obj = x, element = element, elements = elements, 
        p.adjust.method = p.adjust.method, H1.test = H1.test)
    thresh <- NULL
    col.range <- c(-0.5, 0.5)
    if (H1.test != "greater") {
        thresh <- c(thresh, -abs(log10(p.threshold)))
        col.range[1] <- -5
    }
    if (H1.test != "lower") {
        thresh <- c(thresh, abs(log10(p.threshold)))
        col.range[2] <- 5
    }
    cols <- .make.col.obj(n = ncolors, min.col = "red", max.col = "blue", 
        threshold = thresh, range = col.range, threshcol = threshcol, colsqueeze=colsqueeze)
    layout(matrix(1:4, nrow = 2), widths = c(0.3, 0.7), heights = c(0.3, 
        0.7))
    .plot.caption(col.obj = cols, main = element, p.adjust.method = p.adjust.method, 
        thresh.lines = thresh)
    if ("phylo" %in% class(tree)) {
        .plot.phylo(tree, species, horizontal = TRUE)
    }
    else {
        frame()
    }
    if ("phylo" %in% class(tree)) {
        .plot.phylo(tree, species)
    }
    else {
        frame()
    }
    .plot.matrix(pmatrix = stats, species = species, elements = elements, 
        col.obj = cols, ...)
    layout(1)
    par(op)
    ans <- list(tree = tree, species = species, elements = elements, 
        stats = stats)
    return(invisible(ans))
}
