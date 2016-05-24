part.data.sp <-
function(data, data.sp, part, p)
{
    mat <- sapply(which(part == -p), l.actu.data.sp, d=data, d.sp=data.sp)
    pos <- as.integer(levels(as.factor(mat[1,])))
    z.actu <- data.sp$z[pos] - as.vector(tapply(mat[2,], mat[1,], sum))
    n.actu <- data.sp$n[pos] - as.vector(tapply(mat[3,], mat[1,], sum))
    data.sp <- actu.data.sp(pos, z.actu, n.actu, data.sp)

    return(data.sp)
}
