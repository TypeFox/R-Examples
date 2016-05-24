


dcomplnorm<-function(x, spec, sigma=1, theta=1, log=FALSE, ...)
{
        f2<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        F2<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        ff<-function (z) {f2(z, ...)}
        f2p<-grad(ff,theta)

	mu<-log(theta) + sigma^2 + theta*sigma^2*f2p / f2(theta, ...)
        f1<-function (z) {dlnorm(z,meanlog=mu,sdlog=sigma)}
        F1<-function (z) {plnorm(z,meanlog=mu,sdlog=sigma)}
        f1L<-function (z) {dlnorm(z,meanlog=mu,sdlog=sigma,log=TRUE)}
        F1L<-function (z) {plnorm(z,meanlog=mu,sdlog=sigma,log.p=TRUE)}
	phi<-( f1(theta)*(1-F2(theta, ...) ) / ( f2(theta, ...)*F1(theta) ) )

	pdf<-x
        pdf[x<=theta&log==FALSE]<-f1(x[x<=theta]) / ((1+phi)*F1(theta))
        pdf[x<=theta&log==TRUE]<-f1L(x[x<=theta])-log(1+phi)-F1L(theta)
        pdf[x>theta&log==FALSE]<-phi*f2(x[x>theta], ...)/((1+phi)*(1-F2(theta, ...)))
        pdf[x>theta&log==TRUE]<-log(phi)+log(f2(x[x>theta], ...))-log(1+phi)-log(1-F2(theta, ...))
	return(pdf)

}



pcomplnorm<-function(x, spec, sigma=1, theta=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        f2<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        F2<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        ff<-function (z) {f2(z, ...)}
        f2p<-grad(ff,theta)

	mu<-log(theta) + sigma^2 + theta*sigma^2*f2p / f2(theta, ...)
        f1<-function (z) {dlnorm(z,meanlog=mu,sdlog=sigma)}
        F1<-function (z) {plnorm(z,meanlog=mu,sdlog=sigma)}
        F1L<-function (z) {plnorm(z,meanlog=mu,sdlog=sigma,log.p=TRUE)}
	phi<-( f1(theta)*(1-F2(theta, ...) ) / ( f2(theta, ...)*F1(theta) ) )

	cdf<-x
        cdf[x<=theta&log.p==FALSE&lower.tail==TRUE]<-F1(x[x<=theta]) / ((1+phi)*F1(theta))
        cdf[x<=theta&log.p==FALSE&lower.tail==FALSE]<-1-F1(x[x<=theta]) / ((1+phi)*F1(theta))
        cdf[x<=theta&log.p==TRUE&lower.tail==TRUE]<-F1L(x[x<=theta])-log(1+phi)-F1L(theta)
        cdf[x<=theta&log.p==TRUE&lower.tail==FALSE]<-log((1+phi)*F1(theta)-F1(x[x<=theta]))-log(1+phi)-F1L(theta)
        cdf[x>theta&log.p==FALSE&lower.tail==TRUE]<-1/(1+phi) + (phi/(1+phi))*(F2(x[x>theta], ...)-F2(theta, ...))/(1-F2(theta, ...))
        cdf[x>theta&log.p==FALSE&lower.tail==FALSE]<-phi/(1+phi) - (phi/(1+phi))*(F2(x[x>theta], ...)-F2(theta, ...))/(1-F2(theta, ...))
        cdf[x>theta&log.p==TRUE&lower.tail==TRUE]<--log(1+phi)+log(1 + phi*(F2(x[x>theta], ...)-F2(theta, ...))/(1-F2(theta, ...)))
        cdf[x>theta&log.p==TRUE&lower.tail==FALSE]<-log(phi)-log(1+phi)+log(1 - (F2(x[x>theta], ...)-F2(theta, ...))/(1-F2(theta, ...)))
	return(cdf)
}



qcomplnorm<-function(p, spec, sigma=1, theta=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        f2<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        F2q<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}
        F2<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        ff<-function (z) {f2(z, ...)}
        f2p<-grad(ff,theta)

	mu<-log(theta) + sigma^2 + theta*sigma^2*f2p / f2(theta, ...)
        f1<-function (z) {dlnorm(z,meanlog=mu,sdlog=sigma)}
        F1q<-function (z) {qlnorm(z,meanlog=mu,sdlog=sigma)}
        F1<-function (z) {plnorm(z,meanlog=mu,sdlog=sigma)}
	phi<-( f1(theta)*(1-F2(theta, ...) ) / ( f2(theta, ...)*F1(theta) ) )

	qf<-p
        qf[p<=1/(1+phi)]<-F1q( p[p<=1/(1+phi)]*(1+phi) * F1(theta) )
        qf[p>1/(1+phi)]<-F2q( (p[p>(1/(1+phi))]*(1+phi)-1)/phi * (1-F2(theta, ...))+ F2(theta, ...), ...)
	return(qf)
}



rcomplnorm<-function(n, spec, sigma=1, theta=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qcomplnorm(u, spec, sigma=sigma, theta=theta, ...)
	return(sf)
}
