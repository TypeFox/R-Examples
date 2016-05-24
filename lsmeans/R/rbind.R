# rbind method for ref.grid objects

rbind.ref.grid = function(..., deparse.level = 1, adjust = "mvt") {
    objs = list(...)
    if (!all(sapply(objs, inherits, "ref.grid")))
        stop("All objects must inherit from 'ref.grid'")
    bhats = lapply(objs, function(o) o@bhat)
    bhat = bhats[[1]]
    if(!all(sapply(bhats, function(b) (length(b) == length(bhat)) 
                   && (sum((b - bhat)^2, na.rm = TRUE) == 0))))
        stop("All objects must have the same fixed effects")
    Vs = lapply(objs, function(o) o@V)
    V = Vs[[1]]
    if(!all(sapply(Vs, function(v) sum((v - V)^2) == 0)))
        stop("All objects must have the same covariances")
    obj = objs[[1]]
    linfcts = lapply(objs, function(o) o@linfct)
    obj@linfct = do.call(rbind, linfcts)
    bnms = unlist(lapply(objs, function(o) o@misc$by.vars))
    grids = lapply(objs, function(o) o@grid)
    gnms = unique(c(bnms, unlist(lapply(grids, names))))
    gnms = setdiff(gnms, c(".wgt.", ".offset.")) # exclude special names
    grid = data.frame(.tmp. = seq_len(n <- nrow(obj@linfct)))
    for (g in gnms)
        grid[[g]] = rep(".", n)
    grid$.tmp. = NULL
    n.before = 0
    for (g in grids) {
        rows = n.before + seq_along(g[[1]])
        n.before = max(rows)
        for (nm in names(g))
            grid[rows, nm] = as.character(g[[nm]])
    }
    avgd.over = unique(unlist(lapply(objs, function(o) o@misc$avgd.over)))
    attr(avgd.over, "qualifier") = " some or all of"
    obj@grid = grid
    update(obj, pri.vars = gnms, by.vars = NULL, adjust = adjust,
           famSize = round((1 + sqrt(1 + 8*n)) / 2, 3),
           estType = "contrast", infer = c(FALSE, TRUE),
           avgd.over = avgd.over)
}


### Subset a reference grid

"[.ref.grid" = function(x, i, adjust = "mvt", ...) {
    x@linfct = x@linfct[i, , drop=FALSE]
    x@grid = x@grid[i, , drop = FALSE]                  
    x = update(x, pri.vars = names(x@grid), 
           famSize = 2, adjust = adjust)
    x@misc$by.vars = NULL
    x
}

