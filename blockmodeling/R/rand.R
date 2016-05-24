"rand" <-
function (tab) #Hubert & Arabie
{
    n <- sum(tab)
	1 + (sum(tab^2) - (sum(rowSums(tab)^2) + sum(colSums(tab)^2))/2)/choose(n, 2)
}

