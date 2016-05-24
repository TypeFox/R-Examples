hFDR.adjust <- function(unadjp, tree.el, alpha = 0.05) {
    # If user does not name unadjusted p-values or tree nodes,
    # assume i^th element of unadjp corresponds to the i^th
    # row / column of tree.
    if(is.null(names(unadjp))) {
        names(unadjp) <- 1:length(unadjp)
    }
    if(!all(names(unadjp) %in% unique(as.vector(tree.el)))) {
        stop("Names of elements in unadjp do not match names of tree nodes")
    }

    p.vals <- data.frame(unadjp, adjp = NA)
    hyp.tree.unadjusted <- new("hypothesesTree", tree = tree.el,
                               p.vals = p.vals, alpha = alpha)

    # Check to see if possible to descend from the root node (if not
    #  significant reject no hypotheses).
    root <- FindRoot(hyp.tree.unadjusted@tree)
    if(hyp.tree.unadjusted@p.vals[root, 'unadjp'] > alpha){
        warning("Root hypothesis p-value equal to ", hyp.tree.unadjusted@p.vals[root, 'unadjp'], ". Fail to reject any hypotheses, terminating procedure.")
        hyp.tree.unadjusted@p.vals[, 'adj.significance'] <- '-'
        return(hyp.tree.unadjusted)
    }
    
    # Perform correction, and format output
    hyp.tree <- hFDR.internal(hyp.tree.unadjusted)
    hyp.tree@p.vals[root, 'adjp'] <- hyp.tree@p.vals[root, 'unadjp']
    hyp.tree@p.vals[, 'adj.significance'] <- SignificanceStars(alpha, hyp.tree@p.vals[, 'adjp'])
    return(hyp.tree)
}
