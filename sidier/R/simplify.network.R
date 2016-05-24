simplify.network<-function(node.names=NA, modules=NA, coordinates=NA, network=NA, shift = 0.5, max.lwd.edge =2, min.lwd.edge =1, max.vertex.size=4, min.vertex.size=2, label.size=1/2.5, bgcol="white", main="")
{
datos<-cbind(node.names,modules,coordinates)
Length<-1-shift
plot(datos[,3:4],type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="",main=main)

NewCoordNodes<-c()

if(shift!=1)
	{
	for (i in (1:length(unique(modules))))
		{
		sub<-subset(datos,datos[,2]==unique(modules)[i])
		if(is.null(dim(sub[,3:4])))
		cent<-sub[,3:4]
		if(is.null(dim(sub[,3:4]))==FALSE)
		cent<-colMeans(sub[,3:4])
#		points(x=cent[1],y=cent[2])
	
		for (j in 1:nrow(sub))
			{
			nodo1<-as.matrix(sub[j,c(3:4)])
			nodo2<-cent
#			segments(x0=nodo1[1],x1=nodo2[1],y0=nodo1[2],y1=nodo2[2],lwd=lwd.edge)
			edge.slope<-(nodo2[2]-nodo1[2])/(nodo2[1]-nodo1[1])
			edge.angle<-atan(edge.slope)#/(pi/180)
			edge.lengthX<-nodo2[1]-nodo1[1]
			edge.lengthY<-nodo2[2]-nodo1[2]
			xNew=nodo2[1]-edge.lengthX*Length
			yNew=nodo2[2]-edge.lengthY*Length
			NewCoordNodes<-c(NewCoordNodes,sub[j,1],xNew,yNew)
			}
		}
	NewCoordNodes<-matrix(NewCoordNodes,ncol=3,byrow=T)
	colnames(NewCoordNodes)<-c("node","NewX","NewY")


#ESTO GENERE LOS LINKS ENTRE LOS NODOS QUE ESTABAN CONECTADOS EN LA RED ORIGINAL
	for(k in 1:(nrow(network)-1))
	for(l in (k+1):nrow(network))
		{
		node1<-NewCoordNodes[which(row.names(network)[k]==NewCoordNodes[,1]),2:3]
		node2<-NewCoordNodes[which(row.names(network)[l]==NewCoordNodes[,1]),2:3]
		if(network[as.numeric(row.names(network)[k]),as.numeric(row.names(network)[l])]!=0)
		lines(x=c(node1[1],node2[1]),y=c(node1[2],node2[2]))#,lwd=3*EdgeWidth[L3],col=colores[L3])	
		}


	points(x=NewCoordNodes[,2],y=NewCoordNodes[,3],col="black",pch=21,bg=bgcol,cex=min.vertex.size)
	text(x=NewCoordNodes[,2],y=NewCoordNodes[,3],col="black", labels=NewCoordNodes[,1],cex=min.vertex.size*label.size)
	}

if(shift==1)
	{

# ESTO PINTA EL NODO EN EL CENTROIDE, CON TAMAÑO PROPORCIONAL AL 'N'
	vertexsizes<-c()
	CENT<-c()
	for (i in (1:length(unique(modules))))
		{
		sub<-subset(datos,datos[,2]==unique(modules)[i])
		cent<-c(unique(modules)[i],colMeans(sub[,3:4]))
		CENT<-c(CENT,cent)
		vertexsizes<-c(vertexsizes,nrow(sub))
		}
	CENT<-matrix(CENT,ncol=3,byrow=T)
	colnames(CENT)<-c("modules","Xmodule","Ymodule")
	if(length(unique(vertexsizes))==1)
		vs<-rep(min.vertex.size,length(unique(modules)))
	if(length(unique(vertexsizes))!=1)
		vs<-((max.vertex.size-min.vertex.size))/(max(vertexsizes)-min(vertexsizes))*(vertexsizes-min(vertexsizes))+min.vertex.size


#ESTO VA A CALULAR LAS CONEXIONES QUE UNEN LOS MODULOS ENTRE SÍ
	conexiones<-matrix(nrow=length(unique(modules)),ncol=length(unique(modules)))
	colnames(conexiones)<-unique(modules)
	row.names(conexiones)<-unique(modules)
	for (i in 1:(length(unique(modules))-1))
	for (j in (i+1):length(unique(modules)))
		{
		sub1<-subset(datos,datos[,2]==i)
		sub2<-subset(datos,datos[,2]==j)
		conexiones[i,j]<-sum(network[match(sub1[,1],rownames(network)),match(sub2[,1],rownames(network))])
		conexiones[j,i]<-conexiones[i,j]
		}

# ESTO GENERA EL VECTOR ANCHURA DE BANDA PARA UNIR LOS MÓDULOS:

	x1<-c(as.matrix(as.dist(conexiones)))
	if(unique(sort(x1))[1]==0)
	minX1<-unique(sort(x1))[2]

	a<-((max.lwd.edge-min.lwd.edge))/(max(x1)-minX1)*(x1-minX1)+min.lwd.edge
	a2<-as.matrix(as.dist(matrix(a,ncol=ncol(conexiones))))
	a2[which(conexiones==0)]<-0


#ESTO GENERA LOS LINKS ENTRE LOS NODOS QUE ESTABAN CONECTADOS EN LA RED ORIGINAL
	NewCoordNodes<-CENT
	for(i in 1:(nrow(conexiones)-1))
	for(j in (i+1):nrow(conexiones))
		{
		node1<-NewCoordNodes[which(row.names(conexiones)[i]==NewCoordNodes[,1]),2:3]
		node2<-NewCoordNodes[which(row.names(conexiones)[j]==NewCoordNodes[,1]),2:3]
		if(conexiones[as.numeric(row.names(conexiones)[i]),as.numeric(row.names(conexiones)[j])]!=0)
		lines(x=c(node1[1],node2[1]),y=c(node1[2],node2[2]),lwd=a2[i,j])#,lwd=3*EdgeWidth[L3],col=colores	[L3])	
		}

	points(x=CENT[,2],y=CENT[,3],col="black",pch=21,bg=unique(bgcol),cex=vs)
	text(x=CENT[,2],y=CENT[,3],col="black", labels=unique(modules),cex=label.size*vs)
	}
}


