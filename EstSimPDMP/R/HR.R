HR <-
function(dat,t,h=NULL,alpha=1/5,bound=Inf){
	if (is.null(h)){
		h<-1/((length(dat))^(alpha))
		}
	b<-dat; c<-.InvYnBis(dat,max(dat)+1)
	sum( ((1/h)*.ker( (t-b)/h )*c)[b<bound] )
}
