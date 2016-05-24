#ncs: # cifras signif q pone en threshold plot
# Error si hay 0 en off-diagnal 
# Ahora no es necesario tener el mismo orden en las matrices para que lo combine bien
# Cambiar nombre a threshold q aparece feo el numero en el nombre. aqui y en otros q lo hagan OK!
#preparar tb para NA
# Poner error en merge si no sale matriz sin ceros
# Manejar NA:  Comprobar que está bien
# mutationSummary que escupa los alineamientos de nt y de indels
#Meterle que pinte las mutaciones
#Ver pq peta el de Moha, q deberia salir 0.75!

zero.thr <-
function(dis,ptPDF=TRUE,ptPDFname="zero_Network.pdf",cex.label=1,cex.vertex=1, bgcol="white",label.col="black",label=colnames(dis),modules=FALSE,moduleCol=NA,modFileName="Modules_summary.txt",ncs=4,na.rm.row.col=FALSE)
{
#require(igraph)
#require(network)

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

## END removing NA ##

## zero ###

if(length(which(dis==0))!=nrow(dis))
{
DIS<-dis
 M<-matrix(0,nrow(dis),nrow(dis))
for (zero1 in 1:(nrow(dis)-1))
for (zero2 in (zero1+1):nrow(dis))
if(dis[zero1,zero2]==0)
{
 for (zero1 in 1:(nrow(dis)-1))
 for (zero2 in (zero1+1):nrow(dis))
 if(dis[zero1,zero2]==0)
 M[zero1,zero2]<-1
}
dis2<-M
}

## zero ###

G<-graph.adjacency(dis2)
A<-as.network.matrix(dis2)

		if(modules==T)
		{
		comuni<-walktrap.community(G)
		tab1<-matrix(nrow=nrow(dis2),ncol=2)
		tab1<-as.data.frame(tab1)
		tab1[,1]<-label
		tab1[,2]<-comuni$membership
#		colo<-colors()[sample(c(1,23,25:152,203:259,361:657),length(unique(tab1[,2])))]
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
vertis<-plot.network(A)
plot.network(A,coord=vertis,vertex.col=as.matrix(bgcol),label=label,usearrows=0,vertex.cex=2.5*cex.vertex,interactive=F, label.pos=5,label.col=label.col,label.cex=0.8*cex.label,main="Distances equal to zero as links")

if(ptPDF==TRUE)
{
pdf(file=ptPDFname)
plot.network(A,coord=vertis,vertex.col=as.matrix(bgcol),label=label,usearrows=0,vertex.cex=2.5*cex.vertex,interactive=F, label.pos=5,label.col=label.col,label.cex=0.8*cex.label,main="Distances equal to zero as links")
dev.off()
#dev.copy2pdf(file=ptPDFname)
}

}
