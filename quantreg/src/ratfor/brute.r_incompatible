# New Version of Brute Force Algorithm for the Powell Estimator
#
# Roger Koenker:  last revision: 12 February 2008 
#
# Problem:  || min{Ax,c} - b ||_tau  = min!
#
# The matrix H is p by m matrix of the n choose p basis indices h on input.
# When called it is assumed (!!) that U = A[H[,1],]^{-1} and xh = U b[H[,1]].
subroutine brutpow(n,p,m,H,A,b,c,x,tau,U,xh,d,jminz,nflag)

integer n,p,m
double precision x(p),A(n,p),b(n),c(n)
double precision U(p,p),d(p),xh(p)
double precision zero, one,tau,pow,minz,z
integer H(p,m),k,findk,jminz,nflag

PARAMETER(zero = 0.0d0, one = 1.d0)

jminz = 1
minz  = pow(n,p,x,A,b,c,tau)
do j = 2,m {
	k = findk(p,H(1,j),H(1,j-1))
	if(k == 0)
		{nflag = 4; return}
	call pivot(n,p,H(1,j-1),H(k,j),H(k,j-1),A,U,d,xh,nflag)
	if(nflag > 0)
		return
	do i = 1,p{
		xh(i) = b(H(i,j))
		}
	call dgemv('N',p,p,one,U,p,xh,1,zero,x,1)
	z = pow(n,p,x,A,b,c,tau)
	if(z < minz) {
		minz  = z
		jminz = j
		}
	}
return
end

#####################################################

integer function findk(p,h,g)
integer p,k,h(p),g(p)
findk = 0
do k = 1,p{
	if(h(k) != g(k)) 
		{findk = k; break}
	}
return
end
	
