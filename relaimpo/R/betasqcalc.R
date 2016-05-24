"betasqcalc" <-
function (ausgabe, covg, p, variances, rank, diff, rela, var.y) 
{
    # Author and copyright holder: Ulrike Groemping

    #This routine is distributed under GPL version 2 or newer.
    #The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    #program that calculates the contribution according to each variables squared standardized coefficient

    #artificially forced to sum to 100% if rela=T, otherwise given in R^2 units

    #squared standardized beta
    wahrbetasq <- solve(covg[2:(p + 1), 2:(p + 1)], covg[2:(p + 
        1), 1])^2 * diag(covg[2:(p + 1), 2:(p + 1)])/var.y

    # normalize
    if (rela) 
        wahrbetasq <- wahrbetasq/sum(wahrbetasq)

    # ranking
    raengebetasq <- p + 1 - rank(wahrbetasq)

    # pairwise differences
    if (diff & p > 2) 
        diffbetasq <- wahrbetasq[nchoosek(p, 2)[1, ]] - wahrbetasq[nchoosek(p, 
            2)[2, ]]
    if (diff & p == 2) 
        diffbetasq <- wahrbetasq[1] - wahrbetasq[2]

    # output results
    slot(ausgabe, "betasq") <- wahrbetasq
    if (rank) 
        slot(ausgabe, "betasq.rank") <- raengebetasq
    if (diff) 
        slot(ausgabe, "betasq.diff") <- diffbetasq
    return(ausgabe)
}

