pintsc.p<-function(traits,control=NA,n.replicates=1000,N.Pearson=15,plot="Results",tails=1)
{
if(N.Pearson<3) stop("N.Pearson must be higher than 2")

# Defining a INT function taking
# covariances instead of traits as input

# cor_X is the correlation matrix and n the number of individuals in traits (without 'NAs')
int.cov <-function(cor_X,n){
	eig_X<-eigen(cor_X, only.values=TRUE)$values
	d <- eig_X
	p <- length (d)
	INT<-sum((d-1)^2)/(p)
	INT.c<-(INT-((p-1)/n))
	pref="PINTsc = "
	pref2="RelPINTsc = "
	pref3="PINTsc.c = "

	names<-matrix(c(pref,pref2,pref3))
	outs<-matrix(c(round(INT, 3),round((INT/(p-1))*100, 3),round(INT.c, 3)))
	row.names(outs)<-names
	colnames(outs)<-""

	outs
	}

####

# To determine "n" for int.cov
  X<-traits
  nas<-length(unique(which(is.na(X),arr.ind=T)[,1]))
  if(nas>0)
  X<-na.exclude(traits)
  if(nas==0)
  X<-traits

####

# Generating the n.replicates correlation matrices
MAT<-list()
for(i in 1:n.replicates)
	{
	n.variables=ncol(traits)
	mat<-as.dist(matrix(nrow=n.variables,ncol=n.variables))
	NE<-(n.variables-1)*(n.variables/2)
	mat[1:length(mat)]<-rPearson(n=NE,N=N.Pearson)
	mat<-as.matrix(mat)
	diag(mat)<-1

	if(i==1)
	values<-int.cov(cor_X=mat,n=nrow(X))
	if(i>1)
	values<-cbind(values,int.cov(cor_X=mat,n=nrow(X)))

	MAT[[i]]<-mat
	}
colnames(values)<-paste("rep",1:n.replicates)

# Estimating values for the real dataset
REAL<-pintsc(traits,control)

# p-value
Pval<-matrix(nrow=nrow(values))
for (vari in 1:nrow(values))
	{
	p1<-(length(which(values[vari,]<REAL[[vari]])))/length(values[vari,])
	p2<-(length(which(values[vari,]>REAL[[vari]])))/length(values[vari,]) # Just bigger than real
	p3<-2*min(p1,p2)


	if (tails==1)
	Pval[vari,]<-p2
	if (tails==2)
	Pval[vari,]<-p3

	if(Pval[vari,]==0) Pval[vari,]<-paste('<',1/n.replicates,sep="")
	}

# Output
out<-list()
mean<-apply(values,1,mean)
rango<-(paste(apply(values,1,min),mean,sep="-"))
Summary<-as.data.frame(cbind(REAL[1:3],mean,Pval))
colnames(Summary)<-c("Real","Simulation mean","p-value")
N<-c(REAL[4],"-","-")
names(N)<-colnames(Summary)
Summary<-rbind(Summary,N)
row.names(Summary)[4]<-"N"

out[[1]]<-MAT
out[[2]]<-values
out[[3]]<-rango
out[[4]]<-mean
out[[5]]<-REAL
out[[6]]<-Pval
out[[6]]<-Summary

names(out)<-c("Simulated.cor","Simulated.int","Simulated.Range","Simulated.Mean","Real.intsc","Summary")

##################

# Plotting the distribution or the results, if any

if(plot=="Pearson.distribution"|plot=="P")
hist(rPearson(n=100000,N=N.Pearson),breaks=20)

if(plot=="Results"|plot=="R")
	{ 
	dev.new(width=10, height=3)
	layout(matrix(1:3,ncol=3))
	for(iplot in 1:3)
		{
		vari<-iplot #variable to plot
		toplot1<-out$Simulated.int[vari,]
		kk<-hist(as.matrix(toplot1),plot=F)
		maxX<-max(kk$breaks,REAL[[vari]])
		limX1<-min(kk$breaks)-0.05*maxX
		limX2<-maxX+0.05*maxX

		Names<-strsplit(row.names(out$Simulated.int)[vari],split=" =")[[1]][1]


		if(substr(Pval[vari,],1,1)=="<")
		MAIN<-paste(Names,"    p ",Pval[vari,],sep="")
		if(substr(Pval[vari,],1,1)!="<")
		MAIN<-paste(Names,"    p = ",Pval[vari,],sep="")
		 hist(as.matrix(toplot1),add=F,freq=F,breaks=28, xlab=Names, main=MAIN,xlim=c(limX1,limX2))
		abline(v=out$Real.int[vari],lty=2,col="red")
		}
	}

#################


out
}


