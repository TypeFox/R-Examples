findParentIndices <- function(i, indices, fam.id, ind.id, parent.id)
{
    indices[(fam.id==fam.id[i]) & (ind.id==parent.id[i])]
}
