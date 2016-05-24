#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

find.h.up <-
function(t, h.val, h.ranges){
        
	loc	<-	max(which(h.ranges<=t))
	if(loc==length(h.ranges)) loc <- loc-1
	h	<-	h.val[loc]
	return(h)
	}

#----------------------------------------------------------------------------------------#


find.h.down <-
function(t, h.val, h.ranges){
        
	loc	<-	max(which(h.ranges<t))
	if(loc==length(h.ranges)) loc <- loc-1
	h	<-	h.val[loc]
	return(h)
	}

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

find.increasingMLE <-
function(x, w=rep(1, length(x)), delta=rep(1, length(x)), plot=FALSE){

	n		<-	length(x)
	u 		<-	order(x)
	x		<-	x[u]
	w		<-	w[u]
	delta	<-	delta[u]  
	cc 		<-	sum(w)-cumsum(w)
	cc		<-	cc[-length(cc)]
	Si		<-	cumsum(cc*diff(x))
	Si		<-	c(0,Si)
	Si		<-	unique(Si)
	delta.1	<-	delta[1:length(Si)]
	delta.1[length(Si)]	<-	0
	loc		<-	c(which(delta.1==1),length(Si))
	Si		<-	Si[loc]-Si[min(loc)]
	axis	<-	0:(length(loc)-1)        
                                                      
	pts		<-	cbind(c(axis,max(axis)),c(Si,0))
	hpts	<-	chull(pts)
	lcm		<-	sort(hpts)
	lcm		<-	lcm[-length(lcm)]

	if(plot==TRUE){
			plot(pts, pch=19, cex=0.5, main="increasing")
			points(pts[lcm,])
			points(pts[lcm,], type="l")
			}

	y		<-	axis[lcm]
	z		<-	Si[lcm]
	slopes	<-	diff(z)/diff(y)
	mle		<-	c(0,1/slopes)
	xx		<-	x[loc]
	ranges	<-	c(0,xx[lcm])
	H		<-	sapply(x, find.H, mle, ranges)
	h		<-	sapply(x, find.h.up, mle, ranges)
	phi		<-	sum(w*H)-sum(delta*log(h))        

	return(list(mle=mle, ranges=ranges, phi=phi, H=H))
	}

#----------------------------------------------------------------------------------------#

find.decreasingMLE <-
function(x, w=rep(1, length(x)), delta=rep(1, length(x)), plot=FALSE){

	n		<-	length(x)
	u		<-	order(x)
	x		<-	x[u]
	w		<-	w[u]
	delta	<-	delta[u]
	cc		<-	sum(w)-cumsum(c(0,w[-length(x)]))
	Si		<-	cumsum(cc*diff(c(0,x)))
	Si		<-	unique(Si)
	delta.1	<-	delta[1:length(Si)]
	loc		<-	which(delta.1==1)
	Si		<-	c(0,Si[loc])
         
	axis	<-  0:length(loc)
	pts		<-	cbind(c(axis,min(axis)),c(Si,max(Si)))
	hpts	<-	chull(pts)
	gcm		<-	sort(hpts)
	gcm		<-	gcm[-length(gcm)]

	if(plot==TRUE){
			plot(pts, pch=19, cex=0.5, main="decreasing")
			points(pts[gcm,])
			points(pts[gcm,],type="l")
			}

	y		<-	axis[gcm]
	z		<-	Si[gcm]
	slopes	<-	diff(z)/diff(y)
	mle		<-	1/slopes
	xx		<-	x[loc]
	ranges	<-	c(0,xx)[gcm]
	H		<-	sapply(x, find.H, mle, ranges)
	h		<-	sapply(x, find.h.down, mle, ranges)
	phi		<-	sum(w*H)-sum(delta*log(h))        

	return(list(mle=mle, ranges=ranges, phi=phi, H=H))
	}

#----------------------------------------------------------------------------------------#

find.modeMLE <-
function(x, w=rep(1,length(x)), mode, delta=rep(1,length(x)), plot=FALSE){

	m0			<-	mode
	n			<-	length(x)
	u			<-	order(x)
	x			<-	x[u]
	w			<-	w[u]
	delta		<-	delta[u] 
	k0			<-	which(x==mode)[1]
        
	if(plot==TRUE){par(mfrow=c(1,2))}
        
	mle1		<-	find.increasingMLE(pmin(x,m0), w, delta, plot=plot)
	mle2    	<-	find.decreasingMLE(x[x>m0]-m0, w[x>m0], delta[x>m0], plot=plot)
	mle			<-	c(mle1$mle, Inf, mle2$mle)
	ranges		<-	c(mle1$ranges, mle2$ranges+m0)
        
	H			<-	sapply(x, find.H, mle, ranges)
	h			<-	c(sapply(x[1:(k0-1)], find.h.up, mle, ranges), length(x), sapply(x[(k0+1):n], find.h.down, mle, ranges))
	delta[k0]	<-	0
	phi			<-	sum(w*H)-sum(delta*log(h))        
        
	return(list(mle=mle, ranges=ranges, phi=phi, H=H, mode=mode))
	}

#----------------------------------------------------------------------------------------#

find.unimodalMLE <-
function(x, w=rep(1,length(x)), delta=rep(1, length(x)), plot=FALSE){
        
	n		<-	length(x)
	u		<-	order(x) 
	x		<-	x[u]	 
	w		<-	w[u]	 
	delta	<-	delta[u] 
	m		<-	max(which(delta>0))

	phi		<-	rep(0, m)                   
	phi[1]	<-	find.decreasingMLE(x,w,delta)$phi
	for(i in 2:(m-1)){phi[i]	<-	find.modeMLE(x, w, x[i], delta)$phi}        
	phi[m]	<-	find.increasingMLE(x,w,delta)$phi
        
	m0		<-	which(phi==min(phi[-1]))[1]
	if(m0==n){mle	<-	find.increasingMLE(x,w,delta, plot=plot)} 
	if(m0<n) {mle	<-	find.modeMLE(x, w, x[m0], delta, plot=plot)}

	return(list(mle=mle$mle, ranges=mle$ranges, phi=mle$phi,H=mle$H, mode=x[m0]))
	}

#----------------------------------------------------------------------------------------#

find.antimodeMLE <-
function(x, w=rep(1,length(x)), antimode, delta=rep(1,length(x)), plot=FALSE) {

	a0		<-	antimode
	n		<-	length(x)
	u		<-	order(x)
	x		<-	x[u]
	w		<-	w[u]
	delta	<-	delta[u] 
	k0		<-	max(which(x<a0))

	if (plot==TRUE){par(mfrow=c(1,2))}
                                                                
	mle1	<-	find.decreasingMLE(mapply(min, x, x[k0]), w, delta, plot=plot)
	mle2	<-	find.increasingMLE(x[(k0+1):n]-x[k0], w[(k0+1):n], delta[(k0+1):n], plot=plot)
	mle		<-	c(mle1$mle, 0, mle2$mle)
	ranges	<-	c(mle1$ranges, mle2$ranges+x[k0])
	if(k0==n-1){ranges	<-	c(ranges, max(x))}
	H		<-	sapply(x, find.H, mle, ranges)
	h		<-	c(sapply(x[1:k0], find.h.down, mle1$mle, mle1$ranges), sapply(x[(k0+1):n], find.h.up, mle2$mle, mle2$ranges))
	delta[length(delta)]	<-	0
	phi		<-	sum(w*H)-sum(delta*log(h))        
       
	return(list(mle=mle, ranges=ranges, phi=phi, H=H, antimode=antimode))
	}

#----------------------------------------------------------------------------------------#


find.ushapedMLE <-
function(x, w=rep(1,length(x)), delta=rep(1, length(x)), plot=FALSE){
        
	n		<-	length(x)
	u		<-	order(x) 
	x		<-	x[u]	 
	w		<-	w[u]	 
	delta	<-	delta[u] 
	m		<-	max(which(delta>0))

	a0		<-	(c(0,x[-m])+x)/2
	phi		<-	rep(0, m)
	phi[1]	<-	find.increasingMLE(x,w,delta)$phi
	for(i in 2:(m-1)){phi[i]	<-	find.antimodeMLE(x, w, a0[i], delta)$phi}        
	phi[m]	<-	find.decreasingMLE(x,w,delta)$phi

	l0		<-	which(phi==min(phi))[1]
	if(l0==1){mle	<-	find.decreasingMLE(x, w, delta, plot=plot)}
	if(l0>1) {mle	<-	find.antimodeMLE(x, w, a0[l0], delta, plot=plot)}

	return(list(mle=mle$mle, ranges=mle$ranges, phi=mle$phi,H=mle$H, antimode=a0[l0]))
	}

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
