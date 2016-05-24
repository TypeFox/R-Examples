getCanonicalLines <- function(oo,p,type=1,para){
	b<-c()
	k<-c()
	for(i in 0:(para$num.tracks-1)){
		gg<-getGrid(oo$x0,oo$y0,p,BAF=i,type=type,para=para)
		reg<-lm(gg[,2]~gg[,1])
		b<-cbind(b,as.numeric(unlist(reg)[1]))
		k<-cbind(k,as.numeric(unlist(reg)[2]))
	}
	cali<-NULL
	cali$b<-b
	cali$k<-k
	cali$oo<-oo
	return(cali)
}