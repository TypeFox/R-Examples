"prattcalc" <- 
function (ausgabe, covg, p, rank, diff, rela, var.y) 
{
    # Author and copyright holder: Ulrike Groemping

    # This routine is distributed under GPL version 2 or newer.
    # The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    # program that calculates the contribution according Pratts method,
    # i.e. standardized beta*correlation = beta*covariance/vary
    wahrpratt <- solve(covg[2:(p + 1), 2:(p + 1)], covg[2:(p + 
        1), 1]) * covg[2:(p + 1), 1]/var.y

    # normalize
    if (rela) 
        wahrpratt <- wahrpratt/sum(wahrpratt)

    # ranking
    raengepratt <- p + 1 - rank(wahrpratt)

    # pairwise differences
    if (diff & p > 2) 
        diffpratt <- wahrpratt[nchoosek(p, 2)[1, ]] - wahrpratt[nchoosek(p, 
            2)[2, ]]
    if (diff & p == 2) 
        diffpratt <- wahrpratt[1] - wahrpratt[2]

    # output results
    slot(ausgabe, "pratt") <- wahrpratt
    if (rank) 
        slot(ausgabe, "pratt.rank") <- raengepratt
    if (diff) 
        slot(ausgabe, "pratt.diff") <- diffpratt
    return(ausgabe)
}

