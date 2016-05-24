`p2colasr` <-
function (Z, nsim = length(Z)) 
{
    if (Z[1] < 0) {
        p1 = rank(Z)[1]/(nsim + 1)
        ptodos = 1 - (rank(Z)/(nsim + 1))
        if (sum(ptodos[ptodos < p1]) == 0) 
            p2 = p1
        else p2 = p1 + max(ptodos[ptodos < p1])
    }
    else {
        p1 = 1 - rank(Z)[1]/(nsim + 1)
        ptodos = rank(Z)/(nsim + 1)
        if (sum(ptodos[ptodos < p1]) == 0) 
            p2 = p1
        else p2 = p1 + max(ptodos[ptodos < p1])
    }
    return(p2)
}

