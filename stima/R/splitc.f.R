splitc.f <-
function(dataRTA,dataNUM,nodemat,minsplit,minbucket,predtree,type,crit="f2")
{
	#nodemat: columns of nodemat indicate nodemembership of the tree so far grown
	#dataNUM: original data (in original order) with factor variables 
	#dataRTA: factor variables are converted into dummy variables
	#predtree=columnnumbers of predictor variables to be used for splitting
	#type=vector indicating if the variable is factor (1) or numeric (2)
	index1<-1:ncol(nodemat)
	index2<-index1[apply(nodemat,2,sum)>=minsplit]
	nodemat2<-matrix(nodemat[,index2],ncol=length(index2))
	
	critn<--1
	if(crit=="F-value") {critn<-2}
	if(crit=="R2change") {critn<-1}
	if(crit=="f2") {critn<-0}

	npt<-length(predtree)
	nnm2<-nrow(nodemat2)
	mnm2<-ncol(nodemat2)
	ni2<-length(index2)
	nRTA<-nrow(dataRTA)
	mRTA<-ncol(dataRTA)
	nNUM<-nrow(dataNUM)
	mNUM<-ncol(dataNUM)
	nTYP<-length(type)
	fvec<-numeric(4)
	fvec<-.Fortran("rs_splitc",fvec=as.double(fvec),
						as.integer(predtree),
						as.integer(nodemat2),
						as.double(index2),
						as.matrix(dataRTA),
						as.matrix(dataNUM),
						as.integer(type),
						as.integer(npt),
                    as.integer(nnm2),as.integer(mnm2),
                    as.integer(ni2),
                    as.integer(nRTA),as.integer(mRTA),
                    as.integer(nNUM),as.integer(mNUM),
                    as.integer(nTYP),as.integer(minbucket),
						as.integer(critn),package="stima")$fvec
	fvec
}
