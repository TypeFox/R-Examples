single.network <-
function(dis,threshold=NA,ptPDF=TRUE,ptPDFname="Network.pdf",bgcol="white",label.col="black",label=colnames(dis),modules=FALSE,moduleCol=NA,modFileName="Modules_summary.txt",na.rm.row.col=FALSE)
{
 #require (igraph)
 #require (network)

if(is.na(threshold)==TRUE) print("ERROR: No threshold value defined")

## BEGIN na.rm
if(length(which(is.na(dis)))!=0 & na.rm.row.col==FALSE) stop("NA values found")
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
## END na.rm

j<-threshold
dis2<-matrix(1,nrow=nrow(dis),ncol=ncol(dis))
row.names(dis2)<-row.names(dis)
lim<-max(dis)*j
fuera<-which(dis>lim)
dis2[fuera]<-0

G<-graph.adjacency(dis2)
A<-as.network.matrix(dis2)

		if(modules==TRUE)
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
		write.table(file=modFileName,tab1,quote=F,row.names=FALSE)
		}

plot.network(A,vertex.col=as.matrix(bgcol),label=label,usearrows=0,vertex.cex=2.5,interactive=F,label.pos=5,label.col=label.col,label.cex=0.8,main=paste("Threshold=",j,sep=" "))

if(ptPDF==TRUE)
{
pdf(file=ptPDFname)
plot.network(A,vertex.col=as.matrix(bgcol),label=label,usearrows=0,vertex.cex=2.5,interactive=F,label.pos=5,label.col=label.col,label.cex=0.8,main=paste("Threshold=",j,sep=" "))
dev.off()
}

#dev.copy2pdf(file=ptPDFname)

}
