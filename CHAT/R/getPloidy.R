getPloidy <- function(p,type=1,sam.dat,oo){
	nn<-2*type
    model<-1
	if(is.na(nn)){stop('Wrong type input!')}
	pl.overall<-weighted.mean(nn*2^(sam.dat[,4]-oo$y),sam.dat[,7],na.rm=TRUE)
	if(model==1){
		pl.tumor<-(pl.overall-nn*(1-p))/p
	}
	else{
		pl.tumor<-(pl.overall-2*(1-p))/p
	}
	return(c(pl.overall,pl.tumor))
}