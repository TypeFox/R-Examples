is.nearCP <- function(seg.x,seg.y,gg,para){
	# test if a segment is near to a canonical point
	Ng<-dim(gg)[1]
	for(i in 1:Ng){
		if(abs(seg.x-gg[i,1])<=para$std.BAF & abs(seg.y-gg[i,2])<=para$std.LRR){
			return(TRUE)
		}
	}
	return(FALSE)
}