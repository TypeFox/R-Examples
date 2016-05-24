wilson.cc=function(a,b,c,d,n,alpha){
		if (a+b+c+d!=n){return(paste("Caution: a+b+c+d is not equal to n."))}
		if (a+b+c+d==n){
		z=qnorm(1-alpha/2)
		wicc_est=(b-c)/n
		if (wicc_est==0){wicc=1}
		if (wicc_est!=0){
		phidachcc=(a*d-b*c)/sqrt((a+b)*(c+d)*(a+c)*(b+d))
		if (phidachcc=="Inf" |phidachcc=="NaN"){phidachcc=0}

		# solving inequality f with roots u3cc and l3 cc
		f=function(teta) abs(teta-(a+c)/n)-1/(2*n)-z*sqrt((teta*(1-teta)/n))>0
		if (a+c==0){l3cc=0
		u3cc=max(uniroot.all(f, lower = 0, upper = 1, n=100000))}
		if (a+c==n){u3cc=1
		l3cc=uniroot.all(f, lower = 0, upper = 1, n=100000)[1]}
		if (a+c!=n& a+c!=0){
		u3cc=max(uniroot.all(f, lower = 0, upper = 1, n=100000))
		l3cc=uniroot.all(f, lower = 0, upper = 1, n=100000)[1]}
		du3cc=u3cc-(a+c)/n
		dl3cc=(a+c)/n-l3cc

		# solving inequality g with roots u2cc and l2cc
		g=function(teta) abs(teta-(a+b)/n)-1/(2*n)-z*sqrt((teta*(1-teta)/n))>0
		if (a+b==0){l2cc=0
		u2cc=max(uniroot.all(g, lower = 0, upper = 1, n=100000))}
		if (a+b==n){u2cc=1
		l2cc=uniroot.all(g, lower = 0, upper = 1, n=100000)[1]}
		if (a+b!=0 & a+b!=n){
		u2cc=max(uniroot.all(g, lower = 0, upper = 1, n=100000))
		l2cc=uniroot.all(g, lower = 0, upper = 1, n=100000)[1]
		du2cc=u2cc-(a+b)/n
		dl2cc=(a+b)/n-l2cc

		deltacc=sqrt(dl2cc^2-2*phidachcc*dl2cc*du3cc+du3cc^2)
		epsicc=sqrt(du2cc^2-2*phidachcc*du2cc*dl3cc+dl3cc^2)
		wi_lowcc=wicc_est-deltacc
		wi_uppcc=wicc_est+epsicc
		c(wi_lowcc,wi_uppcc)
		}
		}
		cint=c(wi_lowcc,wi_uppcc)
		attr(cint, "conf.level") <- 1-alpha
		rval <- list(conf.int = cint, estimate = wicc_est)
		class(rval) <- "htest"
		return(rval)
		}}
