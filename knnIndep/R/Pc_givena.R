Pc_givena <-
function(i,c,a,N){
	if(length(c) != length(a)){
		stop("c and a should have the same length")
	}
	res = rep(-1,length(c))
	pdi = P_di(i,a,N)
	res[pdi == 0] = 0
	res[pdi != 0] = 1-((P_cge_ale(i,c[pdi != 0]+1,a[pdi != 0],N) - P_cge_ale(i,c[pdi != 0]+1,a[pdi != 0]-1,N)) / pdi[pdi != 0])
	res[a > c] = 0
	return(res)
}
