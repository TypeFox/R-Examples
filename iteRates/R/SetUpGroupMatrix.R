SetUpGroupMatrix <-
function(subTax,VecBin){
	pG<-list()
	pG[[1]]<-PosGroups(subTax,VecBin[1])
	i<-2
	while (i<=length(VecBin)){
	temp<-list()
	remainRow<-remainingTax(subTax,pG[[i-1]])
	
		for (j in 1:length(pG[[i-1]][,1])){
			temp[[j]]<-PosGroups(remainRow[j,],VecBin[i])
				}
		pGt<-list(pG[[i-1]],temp)
		cOmbinStep<-mergelist(pGt)
		pG[[i]]<-cOmbinStep
		i<-i+1
		}
	return(list(groups=pG,VecBin=VecBin))
	}

