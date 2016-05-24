zbrent<-function(func,int,tol){

itmax<-100

eps<-3e-8

a<-int[1]

b<-int[2]

fa<-func(a)

fb<-func(b)

if (fa*fb > 0.0) {
print("Root must be bracketed for zbrent")
return(-1)     }

fc<-fb

iter<-1

while (iter < itmax) {
	if (fb*fc > 0) {
		c<-a
		fc<-fa
		d<-b-a
		e<-d   }
	if (abs(fc) < abs(fb)) {
		a<-b
		b<-c
		c<-a
		fa<-fb
		fb<-fc
		fc<-fa         }
	tol1<-2.0*eps*abs(b)+0.5*tol
	xm<-0.5*(c-b)
	if ((abs(xm) <= tol1) || (fb == 0)) return(b)
	if ((abs(e) >= tol1) && (abs(fa) > abs(fb))) {
		s<-fb/fa
		if (a == c)        {
			p<-2.0*xm*s
			q<-1.0-s   } else {
				q<-fa/fc
				r<-fb/fc
				p<-s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
				q<-(q-1.0)*(r-1.0)*(s-1.0)         }
		if (p > 0.0) q<--q
		p<-abs(p)
		if (2.0*p < min((3.0*xm*q-abs(tol1*q)),abs(e*q)))  {
			e<-d
			d<-p/q   } else {
				d<-xm
				e<-d    }   } else {
						d<-xm
						e<-d }
	a<-b
	fa<-fb
	if (abs(d) >tol1) b<-b+d else b<-b + sign(xm)*tol1
	fb<-func(b)  
	iter<-iter+1 }
	print("Zbrent exceeding maximum iterations")
	return(b)
}
