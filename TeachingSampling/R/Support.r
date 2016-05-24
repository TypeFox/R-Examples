Support <- function(N, n, ID=FALSE){
m       <- matrix(0,choose(N,n),n)
sam     <- matrix(0,choose(N,n),n)
for(i in 1:n)
{
	a=0
 	t=i
	for( r in 1:choose(N,n))
	{
	a=a+1
	B <- choose(N-t,n-i)
	if(a > B)        { a=1;		t=t+1 }
	if( t  > N-n+i) { t=m[r,i-1]+1       }
	m[r,i]=t
	sam[r,i]=ID[t]
	}
}
if (ID==FALSE) {return(m)}
else {return(sam)}
}