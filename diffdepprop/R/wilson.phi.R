wilson.phi=function(a,b,c,d,n,alpha){
if (a+b+c+d!=n){return(paste("Caution: a+b+c+d is not equal to n."))}
if (a+b+c+d==n){
		z=qnorm(1-alpha/2)
		wi_estphi=(b-c)/n
		if (wi_estphi==0){wiphi=1}
		if (wi_estphi!=0){
		if (a*d>b*c){phidachphi=max(a*d-b*c-n/2,0)/sqrt((a+b)*(c+d)*(a+c)*(b+d))}
		else{
		phidachphi=(a*d-b*c)/sqrt((a+b)*(c+d)*(a+c)*(b+d))}
		if (phidachphi=="Inf" |phidachphi=="NaN"){phidachphi=0}

		# solving f with roots u3phi and l3phi
		f=function(teta) abs(teta-(a+c)/n)-z*sqrt((teta*(1-teta)/n))
		u3phi=uniroot.all(f, lower = 0, upper = 1)[2]
		l3phi=uniroot.all(f, lower = 0, upper = 1)[1]
		du3phi=u3phi-(a+c)/n
		dl3phi=(a+c)/n-l3phi

		# solving g with roots u2phi and l2phi
		g=function(teta) abs(teta-(a+b)/n)-z*sqrt((teta*(1-teta)/n))
		u2phi=uniroot.all(g, lower = 0, upper = 1)[2]
		l2phi=uniroot.all(g, lower = 0, upper = 1)[1]
		du2phi=u2phi-(a+b)/n
		dl2phi=(a+b)/n-l2phi

		deltaphi=sqrt(dl2phi^2-2*phidachphi*dl2phi*du3phi+du3phi^2)
		epsiphi=sqrt(du2phi^2-2*phidachphi*du2phi*dl3phi+dl3phi^2)
		wi_lowphi=wi_estphi-deltaphi
		wi_uppphi=wi_estphi+epsiphi
		c(wi_lowphi,wi_uppphi)
		}
		cint=c(wi_lowphi,wi_uppphi)
		attr(cint, "conf.level") <- 1-alpha
		rval <- list(conf.int = cint, estimate = wi_estphi)
		class(rval) <- "htest"
		return(rval)
		}}

