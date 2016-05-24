genizicalc <- function (ausgabe, covg, p, rank, diff, rela, var.y) 
{
    ## genizi (1993, Statistica Sinica)
    covg <- cov2cor(covg)
    covx <- covg[2:(p+1),2:(p+1)]
    eigcovx <- eigen(covx, TRUE)
    rootcovx <- eigcovx$vectors%*%diag(sqrt(eigcovx$values))%*%t(eigcovx$vectors)
    wahrgenizi <- rowSums((rootcovx*matrix(solve(rootcovx, covg[2:(p+1),1]),p,p,byrow=TRUE))^2)
    if (rela) 
        wahrgenizi <- wahrgenizi/sum(wahrgenizi)
    raengegenizi <- p + 1 - rank(wahrgenizi)
    if (diff & p > 2) 
        diffgenizi <- wahrgenizi[nchoosek(p, 2)[1, ]] - wahrgenizi[nchoosek(p, 
            2)[2, ]]
    if (diff & p == 2) 
        diffgenizi <- wahrgenizi[1] - wahrgenizi[2]
    slot(ausgabe, "genizi") <- wahrgenizi
    if (rank) 
        slot(ausgabe, "genizi.rank") <- raengegenizi
    if (diff) 
        slot(ausgabe, "genizi.diff") <- diffgenizi
    return(ausgabe)
}
