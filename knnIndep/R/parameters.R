kr <-
function(r,i0,c){
        epsilon = 2*c-2-i0
		res = rep(-1,length(c))
		res[r==1] = 4*epsilon[r==1]+4
		res[r==2] = 6*epsilon[r==2]^2+6*epsilon[r==2]+2
		res[r==3] = 4*epsilon[r==3]^3
		res[r==4] = epsilon[r==4]^2*(epsilon[r==4]-1)^2
		if(any(res == -1)) stop(paste(r,"(r>4 | r < 1 are impossible values)"))
		return(res)
}
parameters <-
function(r,i0,c,N){
	if(length(c) == 0){return(NULL)}
	comb = rep(-1,length(c))
	d = 2*c-2-i0
	comb[d < 0] = -Inf
	kr = kr(r,i0,c)
	comb[kr <= 0] = -Inf
	comb[N-2*c-1 < 0] = -Inf
	comb[N-4*c+i0+r-1 < 0] = -Inf
	if(any(comb == -1)){
		comb[comb == -1] = 2*lchoose(2*c[comb == -1]-2,i0)+lfactorial(i0)+						#region I
				log(kr[comb == -1])+												#region R
				2*(lfactorial(N-2*c[comb == -1]-1) - lfactorial(N-4*c[comb == -1]+i0+r-1))+		#region IIa + IIb
				lfactorial(N-4*c[comb == -1]+i0+r-1) - lfactorial(N-1)				#region III
	}
	return(exp(comb))
}
