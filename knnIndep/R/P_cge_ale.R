P_cge_ale <-
function(i,c,a,N){
	if(length(a) == 0 || length(c) == 0){return(NULL)}
	res = rep(-1,length(c))
	res[c > floor(N/2)] = 0
	res[a < 0] = 0
#	res[a > c] = 0#NaN
	if(any(a == c)){
		res[a == c] = P_cge_ale(i,a[a==c],a[a==c]-1,N) + P_cge_ale(i,a[a==c]+1,a[a==c],N)-P_cge_ale(i,a[a==c]+1,a[a==c]-1,N)+P_ceq(i,a[a==c],N)
	}
	res[N-2*c+1 < 0] = 0
	res[N-4*c+i+3 <0] = 0
	if(any(res == -1,na.rm=TRUE)){
		res[res == -1 & a <= c] = exp(2*lchoose(2*a[res == -1 & a <= c],i)+lfactorial(i)+ 				#region I
				lfactorial(N-2*c[res == -1 & a <= c]+1)-lfactorial(N-4*c[res == -1 & a <= c]+i+3)+	#region II
				lfactorial(N-2*c[res == -1 & a <= c]+1)-lfactorial(N-1))		#region III
	}
	return(res)
}
