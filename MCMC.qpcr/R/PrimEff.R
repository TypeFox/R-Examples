PrimEff <-
function(data,plot=TRUE) {
	prim=data
	prim[,"gene"]=as.factor(as.character(prim[,"gene"]))
	genes=levels(prim[,"gene"])
	
	if (length(genes)>1) {
		graphcols=ceiling(sqrt(length(genes)))
		graphrows=ceiling(length(genes)/graphcols)
		if (plot==TRUE) par(mfrow=c(graphrows,graphcols))
	  }
	gg=c();ii=c();ee=c();eu=c();el=c()
	for (g in genes)
		{
		s=prim[prim[,"gene"]==g,]
		genlm=lm(cq~log(dna,2),s,na.action="na.omit")
		inter=coef(genlm)[1]
		slop=coef(summary(genlm))[2,1]
		sloperr=coef(summary(genlm))[2,2]
		eff=2^(-1/(slop))
		efflo=round(2^(-1/(slop-sloperr)),2)
		effup=round(2^(-1/(slop+sloperr)),2)
		slop=round(slop,2)
		gg=append(gg,g)
		eff=round(eff,2)
		ee=append(ee,eff)
		eu=append(eu,effup)
		el=append(el,efflo)
		ii=append(ii,round(inter,1))
		if (plot==TRUE){
			plot(cq~log(dna,2),s,main=g,xlab=paste("slope =",slop," intercept =",round(inter,1)), sub=paste("E = ",eff," ( ",round(efflo,2)," - ",round(effup,2)," )",sep=""),mgp=c(2.3,1,0))
			abline(genlm)
		}	
	}
	rr=cbind("gene"=gg,"E"=ee,"E.minus.se"=el,"E.plus.se"=eu,"intercept"=ii)
	row.names(rr)=c()
	rr=data.frame(rr)
	rr$E=as.numeric(as.character(rr$E))
	rr$E.minus.se=as.numeric(as.character(rr$E.minus.se))
	rr$E.plus.se=as.numeric(as.character(rr$E.plus.se))
	rr$intercept=as.numeric(as.character(rr$intercept))
	return(rr)
}
