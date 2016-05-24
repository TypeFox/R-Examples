GCD.big=function(x,y,B){
	q=mpfr(0,B)
	devam=1
    a=mpfr(0,B)
    b=mpfr(0,B)
	a[1]=x
	b[1]=y
	if (x<y){
		a[1]=y
		b[1]=x
	}
	i=2
	while (devam==1){
		q[i-1]=floor(a[i-1]/b[i-1])
		b[i]=a[i-1]-q[i-1]*b[i-1]
		a[i]=b[i-1]      
		if (b[i]==0){
			devam=0
		}
		i=i+1
	}
	k=i-2
	g=b[i-2]
	q=q[1:(i-2)]
	result=list(k=k,q=q,g=g)
	return(result)
}

