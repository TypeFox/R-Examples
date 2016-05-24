create.set.multi <-
function (set, msp.groups) 
{
    set.groups <- obtain.groups(set = set, msp.groups = msp.groups)
    if (length(set.groups) == 0) 
        return(set)
    set.new <- set[setdiff(1:length(set), unlist(set.groups))]
    set.new <- c(set.new, names(set.groups))
}
create.sets.multi <-
function (sets, msp.groups) 
{
    lapply(sets, function(set) {
        create.set.multi(set = set, msp.groups = msp.groups)
    })
}
weight.gsets.with.msprot <-
function (gsets, isets.multi, msp.groups) 
{
    lapply(gsets, function(gset) {
        weight.gset.with.msprot(gset = gset, isets.multi = isets.multi, 
            msp.groups = msp.groups)
    })
}
weight.gset.with.msprot <-
function (gset, isets.multi, msp.groups) 
{
    subunit.genes <- unlist(msp.groups)
    if (length(hit <- which(!is.na(match(gset, subunit.genes)))) == 
        0) {
        return(weight.gset.test(gset = gset, isets = isets.multi))
    }
    gset.subunit.genes <- gset[hit]
    gset.protein.genes <- setdiff(gset, gset.subunit.genes)
    if (length(gset.protein.genes) == 0) 
        stop("One gene set have all subunit genes and NO protein genes, please do not use such gene sets.")
    gset.multi <- create.set.multi(set = gset, msp.groups = msp.groups)
    w.protein.genes <- weight.gset.test(gset = gset.multi, isets = isets.multi, 
        glist = gset.protein.genes)
    w.subunit.genes <- weight.gset.test(gset = gset.multi, isets.multi, 
        glist = gset.subunit.genes)
    c(w.protein.genes, w.subunit.genes)
}
obtain.groups <-
function (set, msp.groups) 
{
    groups.return <- lapply(msp.groups, function(groups) {
        if (length(ind <- which(!is.na(match(set, groups)))) > 
            0) 
            ind
        else NULL
    })
    groups.return[which(!unlist(lapply(groups.return, is.null)))]
}
