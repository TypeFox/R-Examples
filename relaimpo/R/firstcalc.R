"firstcalc" <-
function (ausgabe, p, variances, rank, diff, rela, var.y) 
{
    # Author and copyright holder: Ulrike Groemping

    # This routine is distributed under GPL version 2 or newer.
    # The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    # program that calculates the contribution according to each variables marginal contribution (first in the model)

    # artificially forced to sum to 100% if rela=T, otherwise given in R^2 units
    wahr1st <- as.numeric(variances[[1]] - variances[[2]])

    # normalize
    if (rela) 
        wahr1st <- wahr1st/sum(wahr1st)
    else wahr1st <- wahr1st/var.y

    # ranking
    raenge1st <- p + 1 - rank(wahr1st)

    # pairwise differences
    if (diff & p > 2) 
        diff1st <- wahr1st[nchoosek(p, 2)[1, ]] - wahr1st[nchoosek(p, 
            2)[2, ]]
    if (diff & p == 2) 
        diff1st <- wahr1st[1] - wahr1st[2]

    # output results
    slot(ausgabe, "first") <- wahr1st
    if (rank) 
        slot(ausgabe, "first.rank") <- raenge1st
    if (diff) 
        slot(ausgabe, "first.diff") <- diff1st
    return(ausgabe)
}

