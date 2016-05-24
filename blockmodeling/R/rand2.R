"rand2" <-
function (clu1,clu2) #Hubert & Arabie
{
	tab<-table(clu1,clu2)
	1 + (sum(tab^2) - (sum(rowSums(tab)^2) + sum(colSums(tab)^2))/2)/choose(sum(tab), 2)
}

