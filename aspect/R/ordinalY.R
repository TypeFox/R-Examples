`ordinalY` <-
function(d,yy,r,itermax=100,eps=1e-6,verbose=0) {
rr<-min(r,dim(yy))
if (r == 1)	return(singOrd(d,yy,rr,itermax=itermax,eps=eps,verbose=verbose))
	else return(multOrd(d,yy,r,itermax=itermax,eps=eps,verbose=verbose))
}

