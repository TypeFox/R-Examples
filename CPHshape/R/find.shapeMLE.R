find.shapeMLE <-
function(x, delta=rep(1, length(x)), type="increasing", plot=FALSE){
	
	if(type=="increasing"){
		mle		<-	find.increasingMLE(x, w=rep(1, length(x)), delta, plot)
		mode	<-	NA
		} 
	if(type=="decreasing"){
		mle		<-	find.decreasingMLE(x, w=rep(1, length(x)), delta, plot)
		mode	<-	NA
		} 
	if(type=="unimodal"){
		mle		<-	find.unimodalMLE(x, w=rep(1, length(x)), delta, plot)
		mode	<-	mle$mode
		} 
	if(type=="ushaped"){
		mle		<-	find.ushapedMLE(x, w=rep(1, length(x)), delta, plot)
		mode	<-	mle$antimode
		} 
	
	res			<-	list(beta=NA, h.range=mle$ranges, h.val=mle$mle, phi=mle$phi, H=mle$H, mode=mode, type=type)
	class(res)	<-	"CPHshape"
	return(res)		

	
	} # find.shapeMLE 

