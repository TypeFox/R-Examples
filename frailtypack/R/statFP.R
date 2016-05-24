

statFP <- function(data, indices, fit, dataset, LPcond, LPmarg, groupe, stimeboot, statusboot, ties, tau, marginal, cindex) {
	print("Bootstrap running ...")
	dataset$ligne <- 1:nrow(dataset)   
	databoot <- NA      
	for(i in 1:length(indices)) {
		datatmp<-cbind(dataset[groupe==indices[i],],groupeboot=i)
		databoot <- rbind(databoot,datatmp)
	}
	databoot <- databoot[! is.na(databoot$groupeboot),]
	databoot[,fit$Names.cluster] <- databoot$groupeboot
	indices.unit <- databoot$ligne
# 	model <- fit$call 
# 	model$data<- databoot
# 	fit.boot <- eval(model)
# 	istop <- fit.boot$istop
# 	LPcond <- fit.boot$linear.pred
# 	if (! is.null(fit.boot$frailty.pred)) LPmarg <- LPcond - fit.boot$frailty.pred[databoot$groupeboot] 
# 	else                                  LPmarg <- LPcond
	BAW.boot <- cindexes.frailty(LPcond[indices.unit], LPmarg[indices.unit], stimeboot[indices.unit], statusboot[indices.unit], databoot$groupe, ties, tau, marginal,cindex)
	
	if (marginal==0){
		CPE.B.Cboot <- BAW.boot$CPE.B.C
		CPE.W.Cboot <- BAW.boot$CPE.W.C
		CPE.O.Cboot <- BAW.boot$CPE.O.C
		
		Cuno.B.Cboot <- BAW.boot$Cuno.B.C
		Cuno.W.Cboot <- BAW.boot$Cuno.W.C
		Cuno.O.Cboot <- BAW.boot$Cuno.O.C
		
		if (cindex==1){
			cindex.B.Cboot <- BAW.boot$cindex.B.C
			cindex.W.Cboot <- BAW.boot$cindex.W.C
			cindex.O.Cboot <- BAW.boot$cindex.O.C
		}
	}else{
		CPE.B.Mboot <- BAW.boot$CPE.B.M
		CPE.W.Mboot <- BAW.boot$CPE.W.M
		CPE.O.Mboot <- BAW.boot$CPE.O.M    
		
		Cuno.B.Mboot <- BAW.boot$Cuno.B.M
		Cuno.W.Mboot <- BAW.boot$Cuno.W.M
		Cuno.O.Mboot <- BAW.boot$Cuno.O.M
		
		if (cindex==1){
			cindex.B.Mboot <- BAW.boot$cindex.B.M
			cindex.W.Mboot <- BAW.boot$cindex.W.M
			cindex.O.Mboot <- BAW.boot$cindex.O.M
		}
	}
	
	if (marginal==0) res <- cbind(CPE.B.Cboot, CPE.W.Cboot, CPE.O.Cboot, Cuno.B.Cboot, Cuno.W.Cboot, Cuno.O.Cboot)
	else res <- cbind(CPE.B.Mboot, CPE.W.Mboot, CPE.O.Mboot, Cuno.B.Mboot, Cuno.W.Mboot, Cuno.O.Mboot)
	if (marginal==0 & cindex==1) res <- cbind(res, cindex.B.Cboot, cindex.W.Cboot, cindex.O.Cboot)  
	if (marginal==1 & cindex==1) res <- cbind(res, cindex.B.Mboot, cindex.W.Mboot, cindex.O.Mboot)
	return(res)
}
