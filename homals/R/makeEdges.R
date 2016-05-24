`makeEdges` <-
function(a) {
  c<-(rowSums(a^2))/2
  n<-dim(a)[1]; lns<-matrix(0,0,8)
  for (i in 1:(n-1)) {
  	for(j in (i+1):n) {
  	    dd<-a[i,]-a[j,]; dc<-c[i]-c[j]; ss<-sum(dd^2)
  	    if (is.nul(ss)) next()
  	    ee<-dc*dd/ss; ff<-c(-dd[2],dd[1])
  		xlw<--Inf; xup<-Inf
  		for (k in (1:n)[-c(i,j)]) {
  		    dd<-a[i,]-a[k,]; dc<-c[i]-c[k]
  		    mum<-sum(dd*ff); mom<-dc-sum(dd*ee)
  		    if (is.nul(mum) & (mom > 0)) {
  		    	xlw<-Inf; xup<--Inf
  		    	}
  		    if (mum>0) xlw<-max(xlw,mom/mum)
  		    if (mum<0) xup<-min(xup,mom/mum)
  			}
  		if (xlw<xup) lns<-rbind(lns,c(i,j,ee,ff,xlw,xup))
  		}
  	}
  return(lns)
}

