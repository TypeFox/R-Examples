perc.thr <-
function(dis,range=seq(0,1,0.01),ptPDF=TRUE,ptPDFname="PercolatedNetwork.pdf",estimPDF=TRUE, estimPDFname="PercThr Estimation.pdf",estimOutfile=TRUE, estimOutName="PercThresholdEstimation.txt",cex.label=1,cex.vertex=1,appendOutfile=TRUE,plotALL=FALSE, bgcol="white",label.col="black",label=colnames(dis),modules=FALSE,moduleCol=NA,modFileName="Modules_summary.txt",ncs=4,na.rm.row.col=FALSE, merge=FALSE, save.distance=FALSE,save.distance.name="DistanceMatrix_Perc.thr.txt")
{

if(length(which(is.na(dis)))!=0 & na.rm.row.col==FALSE) stop("NA values found in your input matrix.")

## removing NA ##
if(length(which(is.na(dis)))!=0 & na.rm.row.col==TRUE)
{
dis<-as.matrix(dis)
	repeat
	{
	conNA<-c()
	for (i in 1:nrow(dis))
	conNA<-c(conNA,length(which(is.na(dis[i,]))))
	Out<-sort(which(conNA==sort(conNA,decreasing=T)[1]),decreasing=T)[1]
	dis<-dis[-Out,-Out]
	if(nrow(dis)==0) stop ("The algorithm could not find a matrix without NA values")
	if(length(which(is.na(dis)))==0) break
	}
}

salida<-matrix(nrow=length(range),ncol=3)
colnames(salida)<-c("Threshold","<s>","Clusters")
## END removing NA ##

## merging nodes ##

if(merge==TRUE)
dis<-mergeNodes(dis)

## END merging nodes ##

for (j in range)
{

print(paste("Threshold value:",j,"  Range to test: from ",min(range)," to ",max(range),sep=""))

dis2<-matrix(1,nrow=nrow(dis),ncol=ncol(dis))
lim<-max(dis)*j
fuera<-which(dis>lim)
dis2[fuera]<-0

G<-graph.adjacency(dis2)
A<-as.network.matrix(dis2)

Res<-clusters(G)
noGrande<-sort(Res$csize)[-length(sort(Res$csize))]
N<-sum(noGrande)
repes<-unique(noGrande[which(duplicated(noGrande))])
if (Res$no>1)
{
if (length(repes)>0)
{
n<-c()
for (i in 1:length(repes))
n<-c(n,length(which(noGrande==repes[i])))
sum1<-repes^2*n

noUnic<-c()
for (i in 1:length(repes))
noUnic<-c(noUnic,which(noGrande==repes[i]))
unicos<-noGrande[-noUnic]
sum2<-unicos^2

SUM<-sum(sum1)+sum(sum2)
S<-SUM/N
}

if (length(repes)==0)
S<-sum(noGrande^2)/N
}
if (Res$no==1)
S<-1
#print(c(j,S,Res$no))

salida[which(range==j),1]<-j
salida[which(range==j),2]<-S
salida[which(range==j),3]<-Res$no

if(is.null(colnames(dis)))
label<-c(1:ncol(dis))

if(plotALL==TRUE)
	{
pdf(file=paste("Threshold=",j,".pdf",sep=""))
		if(modules==T)
		{
		comuni<-walktrap.community(G)
		tab1<-matrix(nrow=nrow(dis2),ncol=2)
		tab1<-as.data.frame(tab1)
		tab1[,1]<-label
		tab1[,2]<-comuni$membership
		colores<-tab1[,2]
		bgcol<-colores
		colo<-colour.scheme(def=moduleCol,N=length(unique(tab1[,2])))
		if(is.character(moduleCol[1])==T)
		colo<-moduleCol
		tab1[which(tab1[,2]==1),3]<-colo[1]
		if(length(unique(tab1[,2]))>1)
		for(i in 2:length(unique(tab1[,2])))
		tab1[which(tab1[,2]==i),3]<-colo[i]
		colnames(tab1)<-c("Node_label","Module","Node_colour")
		bgcol<-tab1[,3]
		}

	plot.network(A,vertex.col=as.matrix(bgcol),label=label,usearrows=0,vertex.cex=2.5*cex.vertex, interactive=F,label.pos=5,label.col=label.col,label.cex=0.8*cex.label,main=paste("Percolation threshold=",round(j,ncs),sep=" "))
dev.off()
	}
}

print("Preparing outfiles")

sal<-salida[,-3]
if(estimPDF==TRUE)
{
pdf(file=estimPDFname)
plot(sal,type="l")
points(sal)
dev.off()
}

if(estimPDF==FALSE)
{
plot(sal,type="l")
points(sal)
}

if(estimOutfile==TRUE)
write.table(salida,file=estimOutName,append=appendOutfile,row.names=FALSE)

out<-list(c())
out[[1]]<-salida
out[[2]]<-salida[(max(which(salida[,2]>1))+1),1]
names(out)<-c("Summary","Estimated Percolation Threshold")

#
j<-salida[(max(which(salida[,2]>1))+1),1]
dis2<-matrix(1,nrow=nrow(dis),ncol=ncol(dis))
row.names(dis2)<-row.names(dis)
lim<-max(dis)*j
fuera<-which(dis>lim)
dis2[fuera]<-0

G<-graph.adjacency(dis2)
A<-as.network.matrix(dis2)

		if(modules==T)
		{
		comuni<-walktrap.community(G)
		tab1<-matrix(nrow=nrow(dis2),ncol=2)
		tab1<-as.data.frame(tab1)
		tab1[,1]<-label
		tab1[,2]<-comuni$membership
		#colo<-colors()[sample(c(1,23,25:152,203:259,361:657),length(unique(tab1[,2])))]
		colo<-colour.scheme(def=moduleCol,N=length(unique(tab1[,2])))
		if(is.character(moduleCol[1])==T)
		colo<-moduleCol

		tab1[which(tab1[,2]==1),3]<-colo[1]
		if(length(unique(tab1[,2]))>1)
		for(i in 2:length(unique(tab1[,2])))
		tab1[which(tab1[,2]==i),3]<-colo[i]
		colnames(tab1)<-c("Node_label","Module","Node_colour")

		bgcol<-tab1[,3]
		out[[3]]<-tab1
		names(out)<-c("Summary","Estimated Percolation Threshold","Module")
		
		write.table(file=modFileName,tab1,quote=F,row.names=FALSE)
		}
vertis<-plot.network(A)
plot.network(A,coord=vertis,vertex.col=as.matrix(bgcol),label=label,usearrows=0,vertex.cex=2.5*cex.vertex,interactive=F, label.pos=5,label.col=label.col,label.cex=0.8*cex.label,main=paste("Percolation threshold=",round(j,ncs),sep=" "))

#if(is.na(j)) stop("\n\nNo percolation threshold found.\n\n")

## WARNING IF percolation threshold is not found:
	if(is.na(j) & length(which(dis==0))!=nrow(dis))
	warning("\n\nPercolation threshold can not be estimated and some of the off-diagonal elements in your matrix are zero. Your distance matrix seems to provide low resolution. You may:\n\n1.- Redefine populations by meging those showing distance values of 0 before percolation threshold estimation. For that use the 'merge=TRUE' option \n\n2.- Represent your original distance matrix using the 'No Isolated Nodes Allowed' method. For that use the 'network.method=\"NINA\"' option.\n\n3.- Represent your original distance matrix using the 'zero' method. For that use the 'network.method=\"zero\"' option.")

if(ptPDF==TRUE)
{
pdf(file=ptPDFname)
plot.network(A,coord=vertis,vertex.col=as.matrix(bgcol),label=label,usearrows=0,vertex.cex=2.5*cex.vertex,interactive=F, label.pos=5,label.col=label.col,label.cex=0.8*cex.label,main=paste("Percolation threshold=",round(j,ncs),sep=" "))
dev.off()
}
if(save.distance==TRUE) write.table(file=save.distance.name,dis)

out
}
