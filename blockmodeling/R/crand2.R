"crand2" <-
function (clu1,clu2) #Hubert & Arabie
{
    tab<-table(clu1,clu2)
    n <- sum(tab)
    sum.ni2 <- sum(choose(rowSums(tab), 2))
    sum.nj2 <- sum(choose(colSums(tab), 2))
    E<- sum.ni2 * sum.nj2 / choose(n, 2)
    return((sum(choose(tab, 2)) - E)/((sum.ni2 + sum.nj2)/2 - E))
}

