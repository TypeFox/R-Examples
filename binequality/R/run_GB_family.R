run_GB_family <- function(ID, hb, bin_min, bin_max, obs_mean, ID_name, quantiles = seq(0.006,0.996, length.out = 1000), modelsToFit = c('GB2','GG','BETA2','DAGUM','SINGMAD','LNO','WEI','GA','LOGLOG','PARETO2')){

	controlList<-list('distributions'=c(GB2,GG,GB2,GB2,GB2,LNO,WEI,GA,GB2,GB2),'distNames'=c('GB2','GG','BETA2','DAGUM','SINGMAD','LNO','WEI','GA','LOGLOG','PARETO2'),'linkList'=list(c(muLink=log, sigmaLink=identity, nuLink=log, tauLink=log),c(muLink=log, sigmaLink=log, nuLink=identity, tauLink=NULL),c(muLink=log, sigmaLink=identity, nuLink=log, tauLink=log),c(muLink=log, sigmaLink=identity, nuLink=log, tauLink=identity),c(muLink=log, sigmaLink=identity, nuLink=identity, tauLink=log),c(muLink=identity, sigmaLink=log, nuLink=NULL, tauLink=NULL),c(muLink=log, sigmaLink=log, nuLink=NULL, tauLink=NULL),c(muLink=log, sigmaLink=log, nuLink=NULL, tauLink=NULL),c(muLink=log, sigmaLink=identity, nuLink=identity, tauLink=identity),c(muLink=log, sigmaLink=identity, nuLink=identity, tauLink=log)),'qFuncs'=c(qGB2,qGG,qGB2,qGB2,qGB2,qLNO,qWEI,qGA,qGB2,qGB2),'linkListq'=list(c(muLink=exp, sigmaLink=identity, nuLink=exp, tauLink=exp),c(muLink=exp, sigmaLink=exp, nuLink=identity, tauLink=NULL),c(muLink=exp, sigmaLink=NULL, nuLink=exp, tauLink=exp),c(muLink=exp, sigmaLink=identity, nuLink=exp, tauLink=NULL),c(muLink=exp, sigmaLink=identity, nuLink=NULL, tauLink=exp),c(muLink=identity, sigmaLink=exp, nuLink=NULL, tauLink=NULL),c(muLink=exp, sigmaLink=exp, nuLink=NULL, tauLink=NULL),c(muLink=exp, sigmaLink=exp, nuLink=NULL, tauLink=NULL),c(muLink=exp, sigmaLink=identity, nuLink=NULL, tauLink=NULL),c(muLink=exp, sigmaLink=NULL, nuLink=NULL, tauLink=exp)),'paramCon'=list(list(NULL,NULL,NULL,NULL,FALSE,FALSE,FALSE,FALSE),list(NULL,NULL,NULL,NULL,FALSE,FALSE,FALSE,FALSE),list(NULL,1,NULL,NULL,FALSE,TRUE,FALSE,FALSE),list(NULL,NULL,NULL,1,FALSE,FALSE,FALSE,TRUE),list(NULL,NULL,1,NULL,FALSE,FALSE,TRUE,FALSE),list(NULL,NULL,NULL,NULL,FALSE,FALSE,FALSE,FALSE),list(NULL,NULL,NULL,NULL,FALSE,FALSE,FALSE,FALSE),list(NULL,NULL,NULL,NULL,FALSE,FALSE,FALSE,FALSE),list(NULL,NULL,1,1,FALSE,FALSE,TRUE,TRUE),list(NULL,1,1,NULL,FALSE,TRUE,TRUE,FALSE)),'freeParams'=list(c(TRUE,TRUE,TRUE,TRUE),c(TRUE,TRUE,TRUE,FALSE),c(TRUE,FALSE,TRUE,TRUE),c(TRUE,TRUE,TRUE,FALSE),c(TRUE,TRUE,FALSE,TRUE),c(TRUE,TRUE,FALSE,FALSE),c(TRUE,TRUE,FALSE,FALSE),c(TRUE,TRUE,FALSE,FALSE),c(TRUE,TRUE,FALSE,FALSE),c(TRUE,FALSE,FALSE,TRUE)))
	
	distFits<-list()
	loop.var<-which(controlList$distNames %in% modelsToFit)
	tstamp<-as.numeric(Sys.time())
	
	for(i in loop.var){
		distFits[[i]]=fitFunc(ID = ID, hb = hb, bin_min = bin_min, bin_max = bin_max, obs_mean = obs_mean, ID_name = ID_name, distribution=controlList$distributions[[i]], distName=controlList$distNames[i], links=controlList$linkList[[i]], qFunc=controlList$qFuncs[[i]], quantiles=quantiles, linksq=controlList$linkListq[[i]], con=gamlss.control(c.crit=0.1,n.cyc=200, trace=FALSE), saveQuants = FALSE, muStart=controlList$paramCon[[i]][[1]], sigmaStart=controlList$paramCon[[i]][[2]], nuStart=controlList$paramCon[[i]][[3]], tauStart=controlList$paramCon[[i]][[4]], muFix=controlList$paramCon[[i]][[5]], sigmaFix=controlList$paramCon[[i]][[6]], nuFix=controlList$paramCon[[i]][[7]], tauFix=controlList$paramCon[[i]][[8]], freeParams=controlList$freeParams[[i]], smartStart=TRUE, tstamp=as.numeric(Sys.time()))
	}#end for i
	
	fitComb<-makeFitComb(distFits)
	
	#param filter (undefined means)
	params<-list()
	counter <- 1
	for(j in loop.var){
	  params[[as.character(modelsToFit[counter])]]<-as.data.frame(distFits[[j]]$parameters)
	  counter <- counter + 1
	}
	
	fitComb.filt<-paramFilt(params,fitComb)
	
	bestMod<-modelAvg(fitList = fitComb.filt, ID = ID_name, nonCon=TRUE)
	bestMod.false<-modelAvg(fitList = fitComb.filt, ID = ID_name, nonCon=FALSE)
	
	waicCon<-bestMod$waic
	
	wbicCon<-bestMod$wbic
	
	waicNoCon<-bestMod.false$waic
	
	wbicNoCon<-bestMod.false$wbic
	
	fitComb.filt<-data.frame(fitComb.filt,waicCon,waicNoCon,wbicCon,wbicNoCon)
	
	dat <- data.frame(ID, hb, bin_min, bin_max, obs_mean)
	colnames(dat)[1] <- ID_name
	#LRT
	fitComb.filt<-LRT(dat = dat, fitComb = fitComb.filt, ID = ID_name)
	fitComb<-LRT(dat = dat, fitComb = fitComb, ID = ID_name)
	
	out <- list("fit" = fitComb, "fit.filter" = fitComb.filt, "best_model" = bestMod, "best_model.filter" = bestMod.false)
	
	return(out)
}