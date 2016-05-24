getSobol <- function(d,Index,Coeff,PhiIJ){
	Names = character(0)
	for (nn in 1:d){
		tmp = combn(d,nn)
		for (ii in 1:length(tmp[1,])){
			Names = append(Names,paste(as.character(tmp[,ii]),collapse = ""))
		}
	}
	Values = rep(0,length(Names))
	tmp = Index/Index
	for (nn in 1:d){tmp[nn,] = nn*tmp[nn,]}
	for (mm in 2:length(Index[1,])){
		currentSet = character(0)
		for (nn in 1:d){
			if (!is.nan(tmp[nn,mm])){currentSet = append(currentSet,tmp[nn,mm])}
		}
		currentSet = paste(currentSet,collapse = "")
		
		for (ii in 1:length(Names)){
			if (Names [ii] == currentSet){Values [ii] = Values [ii] + Coeff[mm]^2*PhiIJ[mm]}
		}	
	}
	res = list(Names = Names, Values = Values)
	return(res)
}