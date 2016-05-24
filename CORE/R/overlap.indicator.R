overlap.indicator <-
function(vstart,vend,wstart,wend){
	lw<-length(wstart)
	lv<-length(vstart)
	z<-cbind(c(wstart,vend),c(rep(0,lw),1:lv),c(1:lw,rep(0,lv)))
	z<-z[order(z[,1]),]
	endbefore<-cummax(z[,2])[order(z[,3])][sort(z[,3])!=0]
	z<-cbind(c(vstart,wend),c(1:lv,rep((lv+1),lw)),c(rep(0,lv),1:lw))
	z<-z[order(z[,1]),]
	startafter<-rev(cummin(rev(z[,2])))[order(z[,3])][sort(z[,3])!=0]	
	return(cbind(endbefore+1,startafter-1))
}
