find.H <-
function(t, h.val, h.ranges){
        
	x	<-        h.ranges
	y	<-        h.val
	yy	<-        pmin(y, max(y[y<Inf]))
	xx	<-        pmin(x,t)
	H	<-        sum(yy*diff(xx))
	return(H)
	}

