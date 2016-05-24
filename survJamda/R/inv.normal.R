inv.normal <-
function(i, zstat)
{
	comb = rep(0,length(zstat[[i[1]]]))
	for (k in 1:length(zstat[[i[1]]])){
		for (m in i)
			comb[k]=comb[k]+zstat[[m]][k]
		comb[k]=comb[k]/sqrt(length(i))
	}
	return(comb)
}

