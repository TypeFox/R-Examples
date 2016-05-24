Return.plot <-
function(model,ci=FALSE,alpha=.05)
{
	RP <- c(1.001,2,10,100,1000)
	p.at <- 1-1/RP
	n <- model$n
	p <- 1:n/(n+1)
	emp <- -log(-log(p))
	zF <- distr(p,model=model,type='q')
	plot(emp,sort(model$x.info[,'x']),pch=19,xlim=c(-log(-log(1/1000)),-log(-log(1-1/1000))),xlab='',ylab='',
	ylim=c(distr(1/1000,model=model,type='q'),c(distr(1-1/1000,model=model,type='q'))),axes=FALSE,type='n')
	axis(1,at=-log(-log(p.at)),cex.axis=.8,las=2,labels=round(RP));axis(2,cex.axis=.8,las=2);box()
	title(xlab='Return Period',ylab='Return Level')
	fun.line <- function(x) distr(exp(-exp(-x)),model=model,type='q')
	curve(fun.line,add=TRUE,col='red')
	if(ci)
	{
		fun.low <- function(x) Q.conf.int(exp(-exp(-x)),model=model,alpha=alpha)[1,]
		fun.up <- function(x) Q.conf.int(exp(-exp(-x)),model=model,alpha=alpha)[3,]
		curve(fun.low,add=TRUE,col='red',lty=4)
		curve(fun.up,add=TRUE,col='red',lty=4)
		legend('topleft',paste('Approx. ',100*(1-alpha),'%-CI',sep=''),col='red',lty=4,bty='n')
	}
	points(emp,sort(model$x.info[,'x']),pch=19,cex=.5)
}