getAmpDel <-
function(oo,sam.dat){
	y<-oo$y0
	tmp<-1:dim(sam.dat)[1]
	amp.tag<-which(!tmp%in%oo$list&sam.dat[,4]>y)
	del.tag<-which(!tmp%in%oo$list&sam.dat[,4]<=y)
	amp<-sum(sam.dat[amp.tag,5])/sum(sam.dat[,5])
	del<-sum(sam.dat[del.tag,5])/sum(sam.dat[,5])
	return(c(amp,del))
}
