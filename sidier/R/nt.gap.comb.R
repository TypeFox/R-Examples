nt.gap.comb <-
function(DISTnuc=NA,DISTgap=NA,alpha=seq(0,1,0.1),method="Corrected",saveFile=TRUE,align=NA,silent=FALSE)
{

#sort matrix and check colnames
DISTnuc<-DISTnuc[order(colnames(DISTnuc)),]
DISTnuc<-DISTnuc[,order(colnames(DISTnuc))]

DISTgap<-DISTgap[order(colnames(DISTgap)),]
DISTgap<-DISTgap[,order(colnames(DISTgap))]

cols1<-paste(colnames(DISTnuc),collapse="")
cols2<-paste(colnames(DISTgap),collapse="")

if(cols1!=cols2) stop("Column names of different matrices are not the same!")

CORnuc<-DISTnuc/(max(DISTnuc))
CORgap<-DISTgap/(max(DISTgap))

if(max(DISTnuc)==0) CORnuc<-matrix(0,ncol=ncol(CORnuc),nrow=nrow(CORnuc))
if(max(DISTgap)==0) CORgap<-matrix(0,ncol=ncol(CORgap),nrow=nrow(CORgap))

if(length(alpha)==1)
{#as a matrix if only one alpha value is defined

if(method=="Uncorrected")
{

	alfa<-alpha
	if(alpha=="info")
	{
	EVENTS<-mutationSummary(align,addExtremes=TRUE)[[2]]
	WInd<-EVENTS[,6]/(EVENTS[,3]+EVENTS[,6])
	alfa<-WInd
	}
	DIST<-(1-alfa)*DISTnuc+alfa*DISTgap
	if(saveFile==TRUE)
	write.table(DIST,file=paste	("DistanceMatrixUncorrectedAlpha",alfa,sep="_"))
	OUTuncor<-DIST
}

if(method=="Corrected")
{
	alfa<-alpha
	if(alpha=="info")
	{
	EVENTS<-mutationSummary(align,addExtremes=TRUE)[[2]]
	WInd<-EVENTS[,6]/(EVENTS[,3]+EVENTS[,6])
	alfa<-WInd
	}
	COR<-(1-alfa)*CORnuc+alfa*CORgap
	if(saveFile==TRUE)
	write.table(COR,file=paste("DistanceMatrixCorrectedAlpha",alfa,sep="_"))
	OUTcor<-COR
}
}#as a matrix if only one alpha value is defined



if(length(alpha)>1)
{						#as a list if a alpha is defined
OUTuncor<-list(c())
OUTcor<-list(c())


if(method=="Uncorrected")
{
	for(ite1 in 1:length(alpha))
	{
	alfa<-alpha[ite1]
	DIST<-(1-alfa)*DISTnuc+alfa*DISTgap
	if(saveFile==TRUE)
	write.table(DIST,file=paste	("DistanceMatrixUncorrectedAlpha",alfa,sep="_"))
	OUTuncor[[ite1]]<-DIST
	}
}

if(method=="Corrected")
{
	for(ite1 in 1:length(alpha))
	{
	alfa<-alpha[ite1]
	COR<-(1-alfa)*CORnuc+alfa*CORgap
	if(saveFile==TRUE)
	write.table(COR,file=paste("DistanceMatrixCorrectedAlpha",alfa,sep="_"))
	OUTcor[[ite1]]<-COR
	}
}
}#as a list if a alpha is defined

if(method=="Uncorrected")
OUT<-OUTuncor

if(method=="Corrected")
OUT<-OUTcor

if(silent==FALSE)
print(OUT)
OUT
}
