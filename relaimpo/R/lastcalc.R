"lastcalc" <-
function (ausgabe, p, variances, rank, diff, rela, var.y) 
{
    # Author and copyright holder: Ulrike Groemping

    # This routine is distributed under GPL version 2 or newer.
    # The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    # program that calculates the contribution according to each variables type III SS

    # artificially forced to sum to 100% if rela=T, otherwise given in R^2 units

    wahrIII <- rev(variances[[p]]) - variances[[p + 1]]

    # normalize
    if (rela) 
        wahrIII <- wahrIII/sum(wahrIII)
    else wahrIII <- wahrIII/var.y

    # ranking
    raengeIII <- p + 1 - rank(wahrIII)

    # pairwise differences
    if (diff & p > 2) 
        diffIII <- wahrIII[nchoosek(p, 2)[1, ]] - wahrIII[nchoosek(p, 
            2)[2, ]]
    if (diff & p == 2) 
        diffIII <- wahrIII[1] - wahrIII[2]

    # output results
    slot(ausgabe, "last") <- wahrIII
    if (rank) 
        slot(ausgabe, "last.rank") <- raengeIII
    if (diff) 
        slot(ausgabe, "last.diff") <- diffIII
    return(ausgabe)
}

