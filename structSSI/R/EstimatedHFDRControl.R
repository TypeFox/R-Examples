EstimatedHFDRControl <- function(hyp.tree) {
    alpha <- hyp.tree@alpha
    parent.ix <- rownames(hyp.tree@p.vals) %in% hyp.tree@tree[, 1]
    n.families.tested <- sum(hyp.tree@p.vals[parent.ix, 'adjp'] < alpha, na.rm = T)
    n.tree.discoveries <- sum(hyp.tree@p.vals[, 'adjp'] < alpha, na.rm = T)
    n.tip.discoveries <- sum(hyp.tree@p.vals[!parent.ix, 'adjp'] < alpha, na.rm = T)

    fdr.tree.est <- min(1, (n.tree.discoveries + n.families.tested) /
        (n.tree.discoveries + 1) * alpha)
    fdr.tip.est <-  min(1, (n.tip.discoveries + n.families.tested) /
        (n.tip.discoveries + 1) * alpha)

    return (list(tree = fdr.tree.est,
                 tip = fdr.tip.est,
                 n.families.tested = n.families.tested,
                 n.tree.discoveries = n.tree.discoveries,
                 n.tip.discoveries = n.tip.discoveries))
}
