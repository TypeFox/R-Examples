couple.var.exp.min <-
function(d, prt, ponderer, var.exp.min, cples)
{
    if (!is.matrix(cples))
    {
        couples <- matrix(0, nrow=1, ncol=2)
        couples[1,1] <- cples[1]
        couples[1,2] <- cples[2]
        cples <- couples
    }

    var.exp <- apply(cples, MARGIN=1, calculer.var.exp, data=d, part=prt, ponderer=ponderer)

    if (min(var.exp) > var.exp.min)
    {
        return(list(candidat=FALSE))
    }
    else
    {
        cples <- cples[which(var.exp == min(var.exp)),]
        
        if (is.matrix(cples))
        {
            cples <- cples[1,]
        }

        return(list(candidat=TRUE, vem=min(var.exp), cple=cples))
    }
}
