hFDR.internal <- function(hyp.tree) {
    tree <- graph.edgelist(hyp.tree@tree)
    tree.el.tmp <- data.frame('parent' = hyp.tree@tree[, 1],
                              'child' = hyp.tree@tree[, 2],
                              stringsAsFactors = F)
    root <- FindRoot(tree.el.tmp)
    children <- tree.el.tmp[which(tree.el.tmp$parent == root), 'child']
    children.p.vals <- hyp.tree@p.vals[children, ]
    adjust <- mt.rawp2adjp(children.p.vals$unadjp, 'BH')
    children.p.vals$adjp <- adjust$adjp[order(adjust$index), 'BH']
    hyp.tree@p.vals[children, 'adjp'] <- children.p.vals$adjp

    rejected <- rownames(children.p.vals)[which(children.p.vals$adjp < hyp.tree@alpha)]

    for(child in rejected){
        subcomp <- subcomponent(tree, child, "out")
        if(length(subcomp) > 1){
            subtree.igraph <- induced.subgraph(graph = tree, vids = subcomp)
            subtree.names <- get.vertex.attribute(subtree.igraph, 'name')
            subtree <- new("hypothesesTree", alpha = hyp.tree@alpha,
                           tree = get.edgelist(subtree.igraph))
            subtree@p.vals <- hyp.tree@p.vals[subtree.names, ]
            hyp.tree@p.vals[subtree.names, ] <- hFDR.internal(subtree)@p.vals
        }
    }
    return(hyp.tree)
}
