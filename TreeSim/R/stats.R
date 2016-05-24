stats<-function(lin,merge1,merge2,statstemp){
	statstemp[1] <- statstemp[1]+abs(lin[merge1]-lin[merge2])         #Colless
	statstemp[2]<-statstemp[2]+log(lin[merge1]+lin[merge2]-1)          #s
	statstemp[3]<-statstemp[3]+(lin[merge1]+lin[merge2])               #Sackin (number leaves below each node equiv to number interior vertices above each node)
	if (lin[merge1]+lin[merge2]==2){
		statstemp[4]<-statstemp[4]+1                                   #cherries
	}
	statstemp
	}