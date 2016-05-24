containment.indicator <-
function(vstart,vend,wstart,wend){
	lw<-length(wstart)
	lv<-length(vstart)
	z<-cbind(c(vend,wend),c(1:lv,rep(0,lw)),c(rep(0,lv),1:lw))
	z<-z[order(z[,1]),]
	endbeforeend<-cummax(z[,2])[order(z[,3])][sort(z[,3])!=0]
	z<-cbind(c(wstart,vstart),c(rep((lv+1),lw),1:lv),c(1:lw,rep(0,lv)))
	z<-z[order(z[,1]),]
	startafterstart<-rev(cummin(rev(z[,2])))[order(z[,3])][sort(z[,3])!=0]	
	return(cbind(startafterstart,endbeforeend))
}
