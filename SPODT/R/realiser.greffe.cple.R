realiser.greffe.cple <-
function(mat.adj, bord, cple, part)
{
    mat.adj <- actu.matrice.adj(mat.adj, cple)

    res <- segments.greffe(cple, bord)

    sgmts.grf <- unlist(res$grf)
    bord <- res$bord
    cl.grf <- c(cple[1], cple[2], bord$id[nrow(bord)])

    part[which(part==as.integer(cl.grf[1])  |  part==as.integer(cl.grf[2]))] <- as.integer(cl.grf[3])

    return(list(part=part, adj=mat.adj, sgmts.grf=sgmts.grf, cl.grf=cl.grf, bord=bord))
}
