carcalc <- function (ausgabe, covg, p, rank, diff, rela, var.y) 
{
    ## code adapted from function carscore of package care
    ## by Verena Zuber and Korbinian Strimmer
    ## needs package corpcor
    corg <- cov2cor(covg)
    Rxy <- corg[2:(p+1),1]
    Rxx.invroot <- mpower(corg[2:(p+1),2:(p+1)], -1/2)
    wahrcar <- as.vector(Rxx.invroot %*% Rxy)^2
    if (rela) 
        wahrcar <- wahrcar/sum(wahrcar)
    raengecar <- p + 1 - rank(wahrcar)
    if (diff & p > 2) 
        diffcar <- wahrcar[nchoosek(p, 2)[1, ]] - wahrcar[nchoosek(p, 
            2)[2, ]]
    if (diff & p == 2) 
        diffcar <- wahrcar[1] - wahrcar[2]
    slot(ausgabe, "car") <- wahrcar
    if (rank) 
        slot(ausgabe, "car.rank") <- raengecar
    if (diff) 
        slot(ausgabe, "car.diff") <- diffcar
    return(ausgabe)
}