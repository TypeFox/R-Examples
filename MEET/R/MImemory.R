MImemory<-function(iicc,training.set){

		require("Matrix")

		Prob 	<-as.numeric(iicc$background)
		q 	<- iicc$q
		D	<- iicc$D
		prob_parella<-iicc$probparella
		correction<-iicc$correction_1rOrdre
		correction1rOrdreval<-iicc$correction_1rOrdre_val
		Mperfil<-iicc$Mperfil
		HXmax<-iicc$Hmax
		H<-iicc$HX

		llindar<-slot(correction1rOrdreval,"sderror")[nrow(training.set)]

		matrixSel<-as(iicc$D>llindar,"sparseMatrix")
		interA<-x<-summary(matrixSel)[[1]]
		interB<-y<-summary(matrixSel)[[2]]

		nucleotids<-c("A","T","C","G")
		
		ErrorHX<-slot(correction1rOrdreval,"sderror")[nrow(training.set)+1]
		VarD<-((4-1)^2)/(2*((nrow(training.set)+1)^2)*(log(2,base=exp(1)))^2)
		ErrorMI<-(2*VarD)^(1/2)
	
		Divergence<-lapply(seq(1, length(interA), 1), function(k){
						mi<-matrix(0,length(nucleotids),length(nucleotids))	
						if (interA[k]==interB[k]){
						   diag(mi)<-sapply(c(1:length(nucleotids)),function(i){
											w<-rbind(as.matrix(training.set[,interA[k]]),nucleotids[i])
											ww<-rbind(as.matrix(training.set[,interB[k]]),nucleotids[i])
											training.set.mes.val<-cbind(w,ww)
											pmX<-probability(training.set.mes.val, Prob)
											pmXY<-joint.probability(training.set.mes.val,Prob,prob_parella)
											H<-switch(iicc$classentropy, "Shannon"=entropy.Shannon(pmX),"Renyi"=entropy.Renyi(pmX,q))
											HXY<-entropy.joint(pmXY,q,iicc)
											D<-switch(iicc$classentropy, "Shannon"=divergence.Shannon(training.set.mes.val,H,HXY,correction1rOrdreval),"Renyi"=divergence.Renyi(training.set.mes.val,pmX,pmXY,q,correction1rOrdreval))
											zzmi<-D[1,1]
										})				
						}else{
							
						mi<-sapply(c(1:length(nucleotids)),function(i){
								   zmi<-sapply(c(1:length(nucleotids)),function(j){							
											w<-rbind(as.matrix(training.set[,interA[k]]),nucleotids[i])
											ww<-rbind(as.matrix(training.set[,interB[k]]),nucleotids[j])
											training.set.mes.val<-cbind(w,ww)
											pmX<-probability(training.set.mes.val, Prob)
											pmXY<-joint.probability(training.set.mes.val,Prob,prob_parella)
											H<-switch(iicc$classentropy, "Shannon"=entropy.Shannon(pmX),"Renyi"=entropy.Renyi(pmX,q))
											HXY<-entropy.joint(pmXY,q,iicc)
											D<-switch(iicc$classentropy, "Shannon"=divergence.Shannon(training.set.mes.val,H,HXY,correction1rOrdreval),"Renyi"=divergence.Renyi(training.set.mes.val,pmX,pmXY,q,correction1rOrdreval))
											zzmi<-D[1,2]
											})
								   zmi
								   })
						   }
						   mi
					})
		list(Divergence=Divergence,interA=interB,interB=interA)				
}							
										
										
										
										
										
										
										
										
										
											
			
				
