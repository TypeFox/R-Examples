ApproxC=function(s,p,depth=3){
	A=cumsum(p[s:1])
	B=p
	C=0
	for(i in 1:depth){
	C=C+sum(B[1:(s-i+1)]*A[i:s])/i
	B=convolve(B,rev(p),type="o")
	}
	1/C
}
#fApproxCold=function(s,p){ C=sum(cumsum(p[s:1])*p)+0.5*sum(convolve(p,rev(p),type="o")[1:(s-1)]*cumsum(p[s:1])[1:(s-1)]); 1/C   }
