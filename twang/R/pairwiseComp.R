pairwiseComp <- function(gbm1 = NULL, wtVec = NULL, i = NULL, data, treat.var, vars, sampW = rep(1, nrow(data)), estimand, compFcn, summaryFcn, treatATT = NULL, na.action="level", onlyES = FALSE, onlyKS=FALSE, summarize = TRUE){

	if(onlyES & onlyKS) stop("Only one of \'onlyES\' and \'onlyKS\' can be TRUE.")

	nTreat <- length(levels(data[,treat.var]))
	if(is.null(wtVec)) pScores <- predict(gbm1, type = "response", n.trees = i)[,,1]
	else pScores <- matrix(1/nTreat, ncol = nTreat, nrow = nrow(data))
	pMat <- pScores
	treatMat <- NULL
	treatLev <- levels(data[,treat.var])
	nTreat <- length(treatLev)
	for(i in 1:ncol(pScores)) treatMat <- cbind(treatMat, data[,treat.var] == treatLev[i])
	
	if(estimand == "ATE"){
		pMat <- 1/pMat
		pMat <- pMat * treatMat		
		wtVec <- rowSums(pMat)
	}
	
	if(estimand == "ATT"){
		tATT <- which.max(treatLev == treatATT)
		hldPMat <- pMat
		for(i in 1:ncol(pMat)) pMat[,i] <- hldPMat[,tATT]/hldPMat[,i]
		pMat[,tATT] <- treatMat[,tATT]
		pMat <- pMat * treatMat
		wtVec <- rowSums(pMat)
	}
	
	nComp <- ifelse(estimand == "ATE", choose(nTreat,2), nTreat-1)
	
	hldKS <- NULL
	
	hldSE <- NULL
	
	if(estimand == "ATT"){	
		for(j in 1:(nComp+1)){
			if(j != tATT){
				subDat <- subset(data, data[,treat.var] %in% treatLev[c(j,tATT)])
				subW <- subset(wtVec, data[,treat.var] %in% treatLev[c(j,tATT)])
				subSampW <- subset(sampW, data[,treat.var] %in% treatLev[c(j,tATT)])
				subTreat <- subset(data[,treat.var], data[,treat.var] %in% treatLev[c(j,tATT)])
				subTreat <- subTreat == treatLev[tATT]
				subDat <- cbind(subTreat, subDat)
				#print(bal.stat(subDat, vars = vars, treat.var="subTreat", w.all = subW, 
				#	sampw = subSampW, get.means = TRUE, na.action = na.action, estimand = estimand, 
				#	multinom = TRUE))
				hldKS <- cbind(hldKS, bal.stat(subDat, vars = vars, treat.var="subTreat", w.all = subW, 
					sampw = subSampW, get.means = TRUE, na.action = na.action, estimand = estimand, 
					multinom = TRUE)$results$ks)		
				hldSE <- cbind(hldSE, bal.stat(subDat, vars = vars, treat.var="subTreat", w.all = subW, 
					sampw = subSampW, get.means = TRUE, na.action = na.action, estimand = estimand, 
					multinom = TRUE)$results$std.eff.sz)			
			}
		}
	}
	else{  ## ie, estimand == "ATE"
		for(j in 2:nTreat){
			for(k in 1:(j-1)){
				subDat <- subset(data, data[,treat.var] %in% treatLev[c(j,k)])
				subW <- subset(wtVec, data[,treat.var] %in% treatLev[c(j,k)])
				subSampW <- subset(sampW, data[,treat.var] %in% treatLev[c(j,k)])
				subTreat <- subset(data[,treat.var], data[,treat.var] %in% treatLev[c(j,k)])
				subTreat <- subTreat == treatLev[j]
				subDat <- cbind(subTreat, subDat)
				#print(bal.stat(subDat, vars = vars, treat.var="subTreat", w.all = subW, 
				#	sampw = subSampW, get.means = TRUE, na.action = na.action, estimand = estimand, 
				#	multinom = TRUE))
				hldKS <- cbind(hldKS, bal.stat(subDat, vars = vars, treat.var="subTreat", w.all = subW, 
					sampw = subSampW, get.means = TRUE, na.action = na.action, estimand = estimand, 
					multinom = TRUE)$results$ks)		
				hldSE <- cbind(hldSE, bal.stat(subDat, vars = vars, treat.var="subTreat", w.all = subW, 
					sampw = subSampW, get.means = TRUE, na.action = na.action, estimand = estimand, 
					multinom = TRUE)$results$std.eff.sz)			
					
			}
		}
	}
	
	hldKS <- apply(abs(hldKS), 2, match.fun(compFcn))
	hldSE <- apply(abs(hldSE), 2, match.fun(compFcn))
	if(summarize){
		hldKS <- apply(matrix(hldKS, nrow=1), 1, match.fun(summaryFcn))
		hldSE <- apply(matrix(hldSE, nrow=1), 1, match.fun(summaryFcn))
	}
	
	if(onlyES) return(hldSE)
	else if(onlyKS) return(hldKS)
	else return(list(ks=hldKS, es = hldSE))
}