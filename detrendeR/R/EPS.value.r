EPS.value = function (rwl, stc = c(5, 2, 1)) 
{
    if (sum(stc) != 8) 
        stop("Site-Tree-Core mask does not sum to 8")
    tree.mask <- 8 - stc[3]
    colnames(rwl) <- substr(colnames(rwl), 1, tree.mask)
    tree.names <- colnames(rwl)
    n.cores <- dim(rwl)[2]
    n.trees <- length(unique(colnames(rwl)))
    n.cores.tree = data.frame(table(tree.names))$Freq
    r.mat = cor(rwl)
    n.tot = 0.5 * n.cores * (n.cores - 1)
    rbar.tot = mean(r.mat[upper.tri(r.mat)], na.rm = TRUE)
    within.tree <- matrix(rep(tree.names, n.cores), ncol = n.cores) == 
        matrix(rep(tree.names, n.cores), ncol = n.cores, byrow = T)
    n.wt <- sum(within.tree[upper.tri(within.tree)], na.rm = TRUE)
    {
        if (n.wt == 0)
            rbar.wt = 0
        else {
            within.tree.r <- within.tree * r.mat
            within.tree.r <- replace(within.tree.r, within.tree == 
                FALSE, NA)
            rbar.wt = mean(within.tree.r[upper.tri(within.tree.r)], 
                na.rm = TRUE)
        }
    }
    n.bt = n.tot - n.wt
    rbar.bt = 1/n.bt * (rbar.tot * n.tot - rbar.wt * n.wt)
    c.eff = (1/n.trees * sum(1/n.cores.tree))^-1
    rbar.eff = rbar.bt/(rbar.wt + (1 - rbar.wt)/c.eff)
    n = n.trees
    eps = (n * rbar.eff)/((n * rbar.eff) + (1 - rbar.eff))
    start = FirstYear(rwl)
    end = LastYear(rwl)
    compos.stats = data.frame(start, end, n.trees, n.cores, n.tot, 
        n.wt, n.bt, rbar.tot, rbar.wt, rbar.bt, c.eff, rbar.eff, 
        eps)
    compos.stats = round(compos.stats, 3)
    return(compos.stats)
}



