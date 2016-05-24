wilson=function(a,b,c,d,n,alpha){
	if (a+b+c+d!=n){return(paste("Caution: a+b+c+d is not equal to n."))}
	if (a+b+c+d==n){
	z=qnorm(1-alpha/2)
		wi_est=(b-c)/n
		if (wi_est==0){wi=1}
		if (wi_est!=0){
		phidach=(a*d-b*c)/sqrt((a+b)*(c+d)*(a+c)*(b+d))
		if (phidach=="Inf" |phidach=="NaN"){phidach=0}


		f=function(teta) abs(teta-(a+c)/n)-z*sqrt((teta*(1-teta)/n))
		u3=uniroot.all(f, lower = 0, upper = 1)[2]
		l3=uniroot.all(f, lower = 0, upper = 1)[1]
		du3=u3-(a+c)/n
		dl3=(a+c)/n-l3

		g=function(teta) abs(teta-(a+b)/n)-z*sqrt((teta*(1-teta)/n))
		u2=uniroot.all(g, lower = 0, upper = 1)[2]
		l2=uniroot.all(g, lower = 0, upper = 1)[1]
		du2=u2-(a+b)/n
		dl2=(a+b)/n-l2

		delta=sqrt(dl2^2-2*phidach*dl2*du3+du3^2)
		epsi=sqrt(du2^2-2*phidach*du2*dl3+dl3^2)
		wi_low=wi_est-delta
		wi_upp=wi_est+epsi
		c(wi_low,wi_upp)
		}
		cint=c(wi_low,wi_upp)
		attr(cint, "conf.level") <- 1-alpha
		rval <- list(conf.int = cint, estimate = wi_est)
		class(rval) <- "htest"
		return(rval)
	}}
