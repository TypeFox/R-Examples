realiser.greffe <-
function(d, bord, part, ponderer, var.exp.min, mat.adj)
{
    if (length(unique(part)) < 3)
    {
        return(list(grf=FALSE))
    }

    cples <- coupler.classes.adj(mat.adj)        
    cvem <- couple.var.exp.min(d, part, ponderer, var.exp.min, cples)

    if (cvem$candidat)
    {
        res <- realiser.greffe.cple(mat.adj, bord, cvem$cple, part)
        
        mat.adj <- unlist(res$adj)
        sgmts.grf <- unlist(res$sgmts.grf)
        cl.grf <- unlist(res$cl.grf)
        bord <- res$bord
        part <- unlist(res$part)

        res <- realiser.greffe(d, bord, part, ponderer, var.exp.min, mat.adj)

        if (res$grf)
        {
            mat.adj <- res$adj
            part <- res$part
            sgmts.grf <- rbind(sgmts.grf, unlist(res$sgmts.grf))
            cl.grf <- rbind(cl.grf, res$cl.grf)
            bord <- res$bord

            return(list(part=part, adj=mat.adj, sgmts.grf=sgmts.grf, cl.grf=cl.grf, bord=bord, grf=TRUE))
        }
        else
        {
            return(list(part=part,adj=mat.adj, sgmts.grf=sgmts.grf, cl.grf=cl.grf, bord=bord, grf=TRUE))
        }
    }
    else
    {
        return(list(grf=FALSE))
    }        
}
