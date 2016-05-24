`plotCuspDensities` <-
function(object, main="Conditional density", ...){
	y  <- object$y #model.response(object$model)
	ab <- object$linear.predictors
	alpha <- ab[,'alpha']
	beta  <- ab[,'beta']
	.negbeta <- beta < 0
	.negalpha <- alpha < 0
	.bifset <- between(alpha, cusp.bifset(beta)[,-1])
	.lowplane <-  .negalpha & !.negbeta & !.bifset
	.higplane <- !.negalpha & !.negbeta & !.bifset
    m = list(.lowplane,.higplane,.negbeta,which(.bifset))
	tmp<-sapply(1:4,function(i){
		d<-if(sum(m[[i]]>0)>2) { density(y[m[[i]]]) } else {d<-density(rep(0,2));d$n=0;d$y[]=0;d$y[1]=.5;d}
		xlim <- range(d$x)
		ylim <- range(d$y)
		plot(d, xlim=xlim, ylim=ylim, main=main, ...);
		rug(y[m[[i]]], ...);
		draw.cusp.bifset(mark=i,...)
	})
}

