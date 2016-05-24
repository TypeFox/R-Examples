remainingTax <-
function(subTax,pG){
	
	#provides the matrix of remaining taxa after removing the taxa from groups pG (a matrix from PosGroups)
	remainMatrix<-array(dim=c(length(pG[,1]),length(subTax)-length(pG[1,])))
	RowsPG<-length(pG[,1])
		for (i in 1:RowsPG){
			remainMatrix[i,]<-subTax[-(which(subTax %in% pG[i,]))]
			}
		return(remainMatrix)			
	}

