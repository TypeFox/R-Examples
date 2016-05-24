aps <- function (dataMatrix, dv, ivlist){

k<-length(ivlist)
numcc<-2**k-1

ivlist <- unlist(ivlist)
ilist<-match(ivlist,colnames(dataMatrix))
ilist<-na.omit(ilist)
ilist<-colnames(dataMatrix)[ilist]
dataMatrix<-na.omit(dataMatrix[,c(dv,ilist)])

## Generate an ID for each independent variable to 2^(k-1). 
ivID <- matrix(nrow=k,ncol=1)
ivID[,1]<-2^(0:(k-1))
rownames(ivID)<-ivlist

## Generate a matrix containing the bit representation of each subset.
PredBitMap<-matrix(0, k, numcc)
x<-proc.time()
for (i in 1:numcc){
  n<-i
  j<-1
  while (n != 0){
    if ((n %% 2) == 1) PredBitMap[j,i]<-1
    n <- n%/%2
    j<-j+1
  }
}
rownames(PredBitMap)<-ivlist

apsBitMap<-matrix(nrow=numcc,ncol=1)
apsBitMap[,1]<-t(order(colSums(PredBitMap[,])))

APSMatrix <- matrix(nrow=numcc,ncol=2)
for (i in 1: numcc){
	formula=paste(dv,"~", sep="") 
	for (j in 1: k){
		bit = PredBitMap[j,i]
		if (bit == 1){
			formula=paste(formula,paste("+",ivlist[[j]], sep=""), sep="")
		} 
	}
      APSMatrix[i,2]<-summary(lm(formula,dataMatrix))$r.squared
	APSMatrix[i,1]<-sum(PredBitMap[,i])
}
APSMatrixc<-APSMatrix[apsBitMap[,1],]

## Use the bitmap matrix to generate row headings.
rnl<-vector(mode="character",length=ncol(PredBitMap))
for (i in 1:ncol(PredBitMap)){
  vars<-as.data.frame(PredBitMap[,i])
  colnames(vars)<-"var"
  vars<-subset(vars,var==1)
  rn<-rownames(vars)[1]
  if (nrow(vars)>1){
    for (k in 2:nrow(vars)){
      rn<-paste(rn,rownames(vars)[k],sep=",") 
    }    
  }
  rnl[i]<-rn
}
 
rownames(APSMatrixc)<-rownames(apsBitMap)<-rnl[apsBitMap[,1]]
colnames(APSMatrixc)<-c("k","R2")
return (list(ivID=ivID, PredBitMap=PredBitMap, apsBitMap=apsBitMap, APSMatrix=APSMatrixc))
}

