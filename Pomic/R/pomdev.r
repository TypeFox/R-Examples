`pomdev` <-
function(object1,object2,eps=10^-30,nrange=1000)
{
	if(!is.numeric(object1)|!is.numeric(object2))
        stop("objects must be numeric ")
 
  if(any(!is.finite(object1))|any(!is.finite(object2)))
        stop("objects must have finite values ")
        
  r1<-range(object1)
	r<-extendrange(r1,f=0.01)
	bw1<-bw.nrd0(object1)

  d2<-density(object2,bw=bw1,from=r[1],to=r[2],n=nrange)
	o2<-d2$y
	w <- o2 < eps
	if (any(w)) o2[w] <- eps
	pdf2 <- approxfun( d2$x, o2, yleft=eps, yright=eps)

  x <- - 2* sum ( log(pdf2(object1)) )
  x
}