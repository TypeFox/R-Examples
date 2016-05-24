 Hypergeometric1F1=function(a, b, z)
 {
	if(any(a<=b, b<=0, z<0)).NotYetImplemented()
	df1=2*b
	df2=2*(a-b)
	# ncp=function(x)(2 *(a* z - b *z + x *b *z))/(x *b)
	# curve(ncp, 0, 100)

	ncp=4*z
	x=a/b-1
	#abline(v=x,h=ncp)
	
	df1.2=df1/2; df2.2=df2/2
	fact=x^(1/2*(-2 + df1))*df1^(df1.2)*df2^(df2.2)*(x*df1 + df2)^(1/2*(-df1 - df2))*exp(-ncp/2) / beta(df1.2, df2.2)
 
	stopifnot(x>=0 && df1>0 && df2>0 && ncp>=0)
	df(x, df1, df2, ncp)/fact
 }
 