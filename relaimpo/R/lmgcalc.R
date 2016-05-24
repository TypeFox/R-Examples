"lmgcalc" <-
function (ausgabe, p, indices, variances, rank, diff, rela, var.y) 
{
    # Author and copyright holder: Ulrike Groemping

    #This routine is distributed under GPL version 2 or newer.
    #The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    #program that calculates the LMG variance decomposition
    #(identical to shapley value)
    wahr <- matrix(0, 1, p)
    for (i in 1:p) {
        #add unconditional variance, i.e. j=0
        wahr[i] = wahr[i] + factorial(p - 1) * variances[[1]]
        for (j in 1:p) {
            # j is index of variances, conditional on j-1 variables
             spalten <- which(indices[[j + 1]] == i, arr.ind = T)[, 
                2]
            # j regressors including xi
            if (j > 1) {
                 # indices of columns with j-1 regressors that do not contain i
                 andere <- setdiff(1:ncol(indices[[j]]), which(indices[[j]] == 
                  i, arr.ind = T)[, 2])
                summe <- sum(variances[[j]][andere]) - sum(variances[[j + 
                  1]][spalten])
            }
            if (j == 1) 
                summe <- -sum(variances[[j + 1]][spalten])
            wahr[i] <- wahr[i] + factorial(j - 1) * factorial(p - 
                j) * summe
        }
    }
    # normalize
    if (rela) 
        wahr <- as.numeric(wahr/sum(wahr))
    else wahr <- as.numeric(wahr/(factorial(p) * var.y))
    # ranks
    raenge <- p + 1 - rank(wahr)
    # pairwise differences
    if (diff & p > 2) 
        difflmg <- wahr[nchoosek(p, 2)[1, ]] - wahr[nchoosek(p, 
            2)[2, ]]
    if (diff & p == 2) 
        difflmg <- wahr[1] - wahr[2]

    # output calculation results
    slot(ausgabe, "lmg") <- wahr
    if (rank) 
        slot(ausgabe, "lmg.rank") <- raenge
    if (diff) 
        slot(ausgabe, "lmg.diff") <- difflmg
    return(ausgabe)
}

