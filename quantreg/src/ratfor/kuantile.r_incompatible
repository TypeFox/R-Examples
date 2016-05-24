# Wrapper to compute several quantiles of a sample of n observations
# Calls  K.C. Kiwiel's version of Floyd and Rivest's select algorithm
# Caveat Emptor!!  The ks need to  be sorted.
subroutine kuantile(k,m,n,x)
integer i,j,k(m),m,n
double precision  x(n)

j = 0
do i = 1,m{
	call dsel05(k(i)-j,n-j,x(j+1))
	j = k(i)
	}
return
end
