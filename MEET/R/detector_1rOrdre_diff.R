detector_1rOrdre_diff <-function(training.set, val.set, iicc)
{
	# missing.fun <- iicc$missing.fun
	Prob 	<- as.numeric(iicc$background)
	q<-iicc$q
	HXmax   <- iicc$HXmax
	Redundancia_corregida<- iicc$Redundancia_corregida
	correction <-iicc$correction_1rOrdre 
	ErrorHX<-iicc$ErrorHX
	Herror<-iicc$Herror
	
	seq.rand <-val.set
	training.set.mes.rand<-rbind(as.matrix(training.set),seq.rand)
	pmX	<- probability(training.set.mes.rand, Prob)
	
	H<-switch(iicc$classentropy, "Shannon"=entropy.Shannon(pmX),"Renyi"=entropy.Renyi(pmX,q))
	
	HX 	<- entropy.corrected( H, ErrorHX, HXmax)

	out_det <- diffInstructions(training.set.mes.rand, HX, HXmax,Herror, Redundancia_corregida)
	
	out<-(1/(out_det+.Machine$double.eps))
}

