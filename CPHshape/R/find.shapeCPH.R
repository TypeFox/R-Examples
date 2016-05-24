find.shapeCPH <-
function(x, z, delta, type="increasing", beta0=rep(1, length(as.matrix(z)[1,])), max.loop=250, eps=1e-5, eps.beta=1e-5, print=FALSE){
	
	if(type=="increasing"){
		mle		<-	find.increasingFULLMLE(x, z, delta, beta0, max.loop, eps, eps.beta, print)
		mode	<-	NA
		} 
	if(type=="decreasing"){
		mle		<-	find.decreasingFULLMLE(x, z, delta, beta0, max.loop, eps, eps.beta, print)
		mode	<-	NA
		} 
	if(type=="unimodal"){
		mle		<-	find.unimodalFULLMLE(x, z, delta, beta0, max.loop, eps, eps.beta, print)
		mode	<-	mle$mode
		} 
	if(type=="ushaped"){
		mle		<-	find.ushapedFULLMLE(x, z, delta, beta0, max.loop, eps, eps.beta, print)
		mode	<-	mle$antimode
		} 
	
	res			<-	list(beta=mle$beta, h.range=mle$h.range, h.val=mle$h.val, phi=mle$phi, H=mle$H, mode=mode, type=type)
	class(res)	<-	"CPHshape"
	return(res)		
		
	} # find.shapeCPH 

