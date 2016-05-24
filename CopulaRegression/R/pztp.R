pztp <-
function(y,lambda){
	dummy<-exp(-lambda)
	#cat(paste("length of dummy: ",length(dummy),"\n"))
	#cat(paste("length of ppois: ",length(ppois(y,lambda)),"\n"))
	#cat(paste("---\n"))
    out<-(ppois(y,lambda)-dummy)/(1-dummy)
    out[y<1]=0
return(out)
}
