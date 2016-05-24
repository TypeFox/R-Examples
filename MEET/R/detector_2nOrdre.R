
detector_2nOrdre <-function (training.set, val.set, iicc) {
	
		Prob 	<- iicc$background
		Prob<-as.numeric(Prob)
		q 	<- iicc$q
		D	<- iicc$D
		prob_parella<-iicc$probparella
			
		correction<-iicc$correction_1rOrdre
		correction1rOrdreval<-iicc$correction_1rOrdre_val
		
		Mperfil<-iicc$Mperfil
		HXmax<-iicc$Hmax
		H<-iicc$HX
		
		training.set.mes.rand<-rbind(as.matrix(training.set),val.set)
		
		pmX     <- probability(training.set.mes.rand, Prob)
		pmXY    <- joint.probability(training.set.mes.rand,Prob,prob_parella)
		HXY	<- entropy.joint(pmXY,q,iicc)
	
		H<-switch(iicc$classentropy, "Shannon"=entropy.Shannon(pmX),"Renyi"=entropy.Renyi(pmX,q))

		Dd<-switch(iicc$classentropy, "Shannon"=divergence.Shannon(training.set.mes.rand,H,HXY,correction1rOrdreval),"Renyi"=divergence.Renyi(training.set.mes.rand,pmX,pmXY,q,correction1rOrdreval))
	
		out<- 1/(sum(abs(D-Dd)*Mperfil)+.Machine$double.eps)
		
}

