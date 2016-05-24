P_cge_aeq <-
function(i,c,a,k,N){
	res = rep(-1,length(c))
	res[N-4*c+i+3 <0] = -Inf
	i0 = i-k
	res[i0 < 0] = -Inf
	kr = kr(k[res == -1],i0[res == -1],a[res == -1])
	res[kr <=0] = -Inf
	res[res == -1] = 2*lchoose(2*a[res == -1]-2,i0[res == -1])+lfactorial(i0[res == -1])+log(kr[res == -1])+lfactorial(N-2*c[res == -1]+1)-lfactorial(N-4*c[res == -1]+i+3)+lfactorial(N-2*c[res == -1]+1) - lfactorial(N-1)
	return(exp(res))
}
