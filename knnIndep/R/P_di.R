P_di <-
function(i,a,N){
	res = rep(-1,length(a))
	if(any(a == floor(N/2))){
		res[a == floor(N/2)] = P_ceq(i,a[a == floor(N/2)],N)
	}
	if(any(res == -1)){
		res[res == -1] = P_cge_ale(i,a[res == -1]+1,a[res == -1],N) - P_cge_ale(i,a[res == -1]+1,a[res == -1]-1,N) + P_ceq(i,a[res == -1],N)
	}
	return(res)
}
