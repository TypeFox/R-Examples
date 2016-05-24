subroutine delcol(r,p,k,z,n,nz)
# the p by p-1 matrix r has been formed from a
# p by p upper-triangular by deleting column k
integer p,k,n,nz; double precision r(p,*), z(n,nz)

integer p1,i,j; double precision a,b,c,s,tau

p1 = p-1
for(i = k; i<p; i=i+1) {
# sweep out r(i+1,i)
	a = r(i,i); b = r(i+1,i)
	if(b==0d0) next
# compute the rotation
	if(dabs(b)>dabs(a))
		{ tau = -a/b; s = 1/dsqrt(1d0+tau*tau); c = s * tau }
	else
		{ tau = -b/a; c = 1/dsqrt(1d0+tau*tau); s = c * tau }

# update r and z
	r(i,i) = c*a - s*b
# just for checking; then r(i+1,i)=0
	r(i+1,i) = s*a + c*b
	for(j = i +1; j<=p1; j = j+1) {
		a = r(i,j); b = r(i+1,j)
		r(i,j) = c*a - s *b
		r(i+1,j) = s*a + c * b
	}
	for(j=1; j<=nz; j = j+1) {
		a = z(i,j); b = z(i+1,j)
		z(i,j) = c*a - s*b
		z(i+1,j) = s*a + c*b
	}
}
return
end
