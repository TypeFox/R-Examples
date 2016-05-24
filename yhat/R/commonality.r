commonality <- function (apsOut){

ivID<-apsOut$ivID
nvar<-length(ivID)
if (nvar < 2)return ("Commonality analysis not conducted. Insufficient number of regressors specified.")
numcc<-2**nvar-1
effectBitMap<-apsOut$PredBitMap
apsBitMap<-apsOut$apsBitMap
commonalityMatrix <- matrix(nrow=numcc,ncol=3)
commonalityMatrix[,1]<-apsOut$APSMatrix[order(apsBitMap),2]

## Use the bitmap matrix to generate the list of R2 values needed.
commonalityList<-vector("list", numcc)
for (i in 1: numcc){
	bit = effectBitMap[1,i]
	if (bit == 1) ilist <-c(0,-ivID[1])
		else  ilist<-ivID[1]
	for (j in 2: nvar){
		bit = effectBitMap[j,i]
		if (bit == 1){	
			alist<-ilist
			blist<-genList(ilist,-ivID[j])
			ilist<-c(alist,blist)
		}
		else ilist<-genList(ilist,ivID[j])
	}
	ilist<-ilist*-1
 	commonalityList[[i]]<-ilist
}

## Use the list of R2 values to compute each commonality coefficient.
for (i in 1: numcc){
	r2list <- unlist(commonalityList[i])
	numlist = length(r2list)
	ccsum=0
	for (j in 1:numlist){
		indexs = r2list[[j]]
		indexu = abs (indexs)
		if (indexu !=0) {
			ccvalue = commonalityMatrix[indexu,1]
			if (indexs < 0)ccvalue = ccvalue*-1
			ccsum=ccsum+ccvalue
		}
	} 
	commonalityMatrix[i,2]=ccsum
}
commonalityMatrix<-commonalityMatrix[,-1]
commonalityMatrix[,1]<-commonalityMatrix[apsBitMap[,1],1]
commonalityMatrix[,2]<-commonalityMatrix[,1]/apsOut$APSMatrix[numcc,2]
colnames(commonalityMatrix)<-c("Coefficient", " % Total")
rownames(commonalityMatrix)<-rownames(apsBitMap)
return (commonalityMatrix)
}

