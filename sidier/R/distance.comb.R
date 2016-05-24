distance.comb <-
function(matrices=NA,alphas=NA,method="Corrected",saveFile=TRUE,na.rm=FALSE)
{
nmat<-length(matrices)
mats<-list()
length(mats)<-nmat
CORmats<-list()
length(CORmats)<-nmat

if(is.na(alphas[1])==FALSE & length(alphas)!=nmat) stop ("The number of matrices (",nmat,") and weights (",length(alphas),") are different!!. Check 'alphas' option.")

if(is.na(alphas[1]))
alphas<-rep(1/nmat,nmat)


print(alphas)
sumNA<-function(a,b,na.rm=TRUE)
	{
	SUMna<-matrix(nrow=nrow(a),ncol=ncol(a))
	colnames(SUMna)<-colnames(a)
	row.names(SUMna)<-row.names(a)
	for (i in 1:nrow(a))
	for (j in 1:nrow(b))
		{
		if(is.na(a[i,j])==F&is.na(b[i,j])==F)
		SUMna[i,j]<-a[i,j]+b[i,j]
		if(na.rm==TRUE)
			{
			if(is.na(a[i,j])==T&is.na(b[i,j])==F)
			SUMna[i,j]<-b[i,j]
			if(is.na(a[i,j])==F&is.na(b[i,j])==T)
			SUMna[i,j]<-a[i,j]
			if(is.na(a[i,j])==T&is.na(b[i,j])==T)
			SUMna[i,j]<-NA
			}
		}
	SUMna
	}

rename<-"no"
for (i in 1:nmat)
	{
	if(is.null(colnames(eval(parse(text=matrices[i])))))
		{
		warning("Undefined row names! Matrices will be combined according to their order in the different matrices.")
		rename<-"yes"
		break
		}
	}

names<-c()
for (i in 1:nmat)
	{
	mat.help<-eval(parse(text=matrices[i]))
		if(rename=="yes")
			{colnames(mat.help)<-c(1:nrow(mat.help)) 
			row.names(mat.help)<-c(1:nrow(mat.help))}
	mat.help<-mat.help[order(colnames(mat.help)),]
	mat.help<-mat.help[,order(colnames(mat.help))]
	mats[[i]]<-mat.help
	CORmats[[i]]<-mats[[i]]/(max(mats[[i]],na.rm=T))
	rm(mat.help)
	names<-c(names,paste(colnames(mats[[i]]),collapse=""))
	}
if(length(names)!=length(which(names[1]==names))) stop("Column names of different matrices are not the same!")

MSUM<-mats[[1]]*alphas[[1]]

for (i in 2:nmat)
	{
	if(method=="Uncorrected")
		{
		MATS<-mats[[i]]*alphas[[i]]
		MSUM<-sumNA(MSUM,MATS,na.rm=na.rm)
		OUTuncor<-MSUM
		row.names(OUTuncor)<-colnames(OUTuncor)
			if(saveFile==TRUE)
			write.table(OUTuncor,file="Weighted_Matrix_Uncorrected")
		}

	if(method=="Corrected")
		{
		MATS<-CORmats[[i]]*alphas[[i]]
		MSUM<-sumNA(MSUM,MATS,na.rm=na.rm)
		OUTcor<-MSUM
		row.names(OUTcor)<-colnames(OUTcor)
			if(saveFile==TRUE)
			write.table(OUTcor,file="Weighted_Matrix_Corrected")
		}
	}

if(method=="Uncorrected")
OUT<-OUTuncor

if(method=="Corrected")
OUT<-OUTcor

print(OUT)

}
