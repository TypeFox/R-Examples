actu.data.sp <-
function(pos, z.actu, n.actu, d.sp)
{
    z.actu <- replace(d.sp$z, pos, z.actu)
    n.actu <- replace(d.sp$n, pos, n.actu)
    d.sp <- data.frame(d.sp$loc, z.actu, n.actu)[which(n.actu != 0),]
    colnames(d.sp) <- c("loc", "z", "n")

    return(d.sp)
}
