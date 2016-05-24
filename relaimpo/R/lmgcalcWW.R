"lmgcalcWW" <-
function (ausgabe, p, indices, variances, rank, diff, rela, var.y, WW, ngroups=NULL) 
{
    # Author and copyright holder: Ulrike Groemping

    #This routine is distributed under GPL version 2 or newer.
    #The text of this license can be found at 
    #http://www.gnu.org/copyleft/gpl.html.

    #program that calculates the LMG variance decomposition
    #(identical to shapley value)
    wahr <- matrix(0, 1, p)
    for (i in 1:p) {
      #add unconditional variance, i.e. j=0
      wahr[i] = wahr[i] + factorWW(p, WW, i, ngroups=ngroups) * variances[[1]]
      for (j in 1:p) {
        # j is index of variances, conditional on j-1 variables
#        spalten <- numeric(0)
        spalten <- NULL
        if(any(indices[[j + 1]] == i)) spalten <- 
                    which(indices[[j + 1]] == i, arr.ind = T)[, 2]
        if (length(spalten) > 0) {
          # j regressors including xi
          if (j > 1) {
            # indices of columns with j-1 regressors that do not contain i
            andere <- 1:ncol(indices[[j]])
            if (any(indices[[j]] == i)) 
                andere <- setdiff(1:ncol(indices[[j]]), 
                               which(indices[[j]] == i, arr.ind = T)[, 2])
            if (length(andere) > 0) {
                if (max(length(spalten),length(andere)) > 1) {
                    # Welche Spalten passen zusammen
                    passen <- matrix(rep(FALSE,length(spalten)*length(andere)), 
                                                             length(spalten))
                    for (r in 1:length(spalten))
                        for (s in 1:length(andere))
                            passen[r,s] <- all(indices[[j]][,andere[s]] %in% 
                                 c(indices[[j+1]][,spalten[r]],i))
                    spalten <- spalten[which(passen == TRUE, arr.ind = T)[,1]]
                    andere  <- andere[which(passen == TRUE, arr.ind = T)[,2]]
                }
                WW.fac  <- numeric(length(spalten))
                for (k in 1:length(spalten)) 
                    WW.fac[k] <- factorWW(p, WW, indices[[j+1]][,spalten[k]], 
                                  indices[[j]][,andere[k]], ngroups=ngroups)
                summe <- sum(variances[[j]][andere] * WW.fac) - 
                              sum(variances[[j+1]][spalten] * WW.fac)
            } else {
                summe <- 0
            }
          }
                if (j == 1) 
                    summe <- -sum(variances[[j + 1]][spalten] * 
                                   factorWW(p, WW, i, ngroups=ngroups))
                wahr[i] <- wahr[i] + summe
            }
        }
    }
    # normalize
    if (rela) {
        wahr <- as.numeric(wahr/sum(wahr))
    } else {
        wahr <- as.numeric(wahr/(factorWW(p, WW, ngroups=ngroups) * var.y))
    }
    # ranks
    raenge <- p + 1 - rank(wahr)
    # pairwise differences
    if (diff & p > 2) 
        difflmg <- wahr[nchoosek(p,2)[1,]] - wahr[nchoosek(p,2)[2,]]
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


