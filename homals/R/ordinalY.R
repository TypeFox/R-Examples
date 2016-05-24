`ordinalY` <-
function(d,y,r,itermax=100,eps=1e-6,verbose=0) {
r<-min(r,dim(y))
if (r == 1)	return(singOrd(d,y,r,itermax=itermax,eps=eps,verbose=verbose))
	else return(multOrd(d,y,r,itermax=itermax,eps=eps,verbose=verbose))
}

