mutation.network <-
function(align=NA,indel.method="MCIC",substitution.model="raw",pairwise.deletion=TRUE,network.method="percolation",range=seq(0,1,0.01),
addExtremes=FALSE,alpha="info",combination.method="Corrected",na.rm.row.col=FALSE, modules=FALSE,moduleCol=NA,modFileName="Modules_summary.txt", save.distance=FALSE, save.distance.name="DistanceMatrix_threshold.txt",silent=FALSE,
bgcol="white", label.col="black",label=NA,label.sub.str=NA,colInd="red", colSust="black",lwd.mut=1,lwd.edge=1.5,cex.mut=1,cex.label=1,cex.vertex = 1,main="",
InScale=1, SuScale=1, legend=NA, legend.bty="o",legend.pos="bottomright")



{
#### ALIGNMENT OF UNIQUE HAPLOTYPES:
#
alignUnique<-GetHaplo(align=align, saveFile =FALSE, format = "fasta", seqsNames = NA,silent=T)
#
#
### BEGIN MUTATION METHODS ###########
#
#SUBSTITUTIONS:
#
SuDist<-as.matrix(dist.dna(x=alignUnique,model=substitution.model,pairwise.deletion=pairwise.deletion))
#
#INDELS
#
#1-SIC
if(indel.method=="SIC")
InDist<-SIC(align=alignUnique, saveFile = F, addExtremes = addExtremes)[[2]]
#
#2-FIFTH
if(indel.method=="FIFTH")
InDist<-FIFTH(align=alignUnique, saveFile = F, addExtremes = addExtremes)
#
#3-BARRIEL
if(indel.method=="BARRIEL")
InDist<-BARRIEL(align=alignUnique, saveFile = F, addExtremes = addExtremes)[[2]]
#
#4-MCIC
if(indel.method=="MCIC")
InDist<-MCIC(align=alignUnique, saveFile = F,silent=T)
#
### END MUTATION METHODS ###########
#
## BEGIN MATRIX COMBINATION
if(sum(as.data.frame(InDist))==0&sum(as.data.frame(SuDist))==0) stop("Incorrect distance matrix. All sequences are identical!")
if(sum(as.data.frame(InDist))==0&sum(as.data.frame(SuDist))!=0) dis<-as.matrix(SuDist)
if(sum(as.data.frame(InDist))!=0&sum(as.data.frame(SuDist))==0) dis<-as.matrix(InDist)
if(sum(as.data.frame(InDist))!=0&sum(as.data.frame(SuDist))!=0)
dis<-nt.gap.comb(DISTgap=InDist, DISTnuc=SuDist, alpha=alpha, method=combination.method, saveFile=F,align=alignUnique,silent=TRUE)
#
## END MATRIX COMBINATION
#
#
## SOME ERRORS
if(length(which(is.na(dis)))!=0 & na.rm.row.col==FALSE) stop("NA values found")
#
#
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

## merging nodes ##

if(length(which(dis==0))!=nrow(dis))
if(merge==TRUE)
dis<-mergeNodes(dis)

## END merging nodes ##
#
#
### BEGIN THRESHOLD ESTIMATION ###

## 1- PERCOLATION THRESHOLD
	if(network.method=="percolation")
	{
	salida<-matrix(nrow=length(range),ncol=2)
	colnames(salida)<-c("Threshold","#Clusters")
		for (j in range)
		{
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

		salida[which(range==j),1]<-j
		salida[which(range==j),2]<-Res$no

		}
	j<-salida[(max(which(salida[,2]>1))+1),1]
	dis2<-matrix(1,nrow=nrow(dis),ncol=ncol(dis))
	row.names(dis2)<-row.names(dis)
	lim<-max(dis)*j
	fuera<-which(dis>lim)
	dis2[fuera]<-0
	}
## 1- END PERCOLATION THRESHOLD
#
#
## 2- NINA THRESHOLD

if(network.method=="NINA")
	{
	salida<-matrix(nrow=length(range),ncol=2)
	colnames(salida)<-c("Threshold","#Clusters")
		for (j in range)
		{

		dis2<-matrix(1,nrow=nrow(dis),ncol=ncol(dis))
		lim<-max(dis)*j
		fuera<-which(dis>lim)
		dis2[fuera]<-0

		G<-graph.adjacency(dis2)
		A<-as.network.matrix(dis2)

		Res<-clusters(G)
		salida[which(range==j),1]<-j
		salida[which(range==j),2]<-Res$no
		}
	j<-salida[min(which(salida[,2]==1)),1]
	dis2<-matrix(1,nrow=nrow(dis),ncol=ncol(dis))
	row.names(dis2)<-row.names(dis)
	lim<-max(dis)*j
	fuera<-which(dis>lim)
	dis2[fuera]<-0
	}

## 2- END NINA THRESHOLD
#
#
## 3- ZERO THRESHOLD
if(network.method=="zero")
	{
	if(length(which(dis==0))==nrow(dis)) stop ("No offdiagonal zeros in your input matrix. Use another nerwork method.")
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
		j<-0
		}
	}
## 3- END ZERO THRESHOLD

if(network.method=="zero")
	{
	DIS<-dis
	for (zero1 in 1:(nrow(dis)-1))
	for (zero2 in (zero1+1):nrow(dis))
	if(dis[zero1,zero2]==0)
		{
		M<-matrix(0,nrow(dis),nrow(dis))
		for (zero1 in 1:(nrow(dis)-1))
		for (zero2 in (zero1+1):nrow(dis))
		if(dis[zero1,zero2]==0)
		M[zero1,zero2]<-1
		}
	dis2<-M
	}


### END ZERO THRESHOLD ###
#
#
### END THRESHOLD ESTIMATION ###
#
#
#write.table(dis2,file="kk_control_combined_dis2.txt")
#
## WARNING IF percolation threshold is not found:
	if(is.na(j) & length(which(dis==0))!=nrow(dis))
	warning("\n\nPercolation threshold can not be estimated and some of the off-diagonal elements in your matrix are zero. Your distance matrix seems to provide low resolution. You may:\n\n1.- Redefine populations by meging those showing distance values of 0 before percolation threshold estimation. For that use the 'merge=TRUE' option \n\n2.- Represent your original distance matrix using the 'No Isolated Nodes Allowed' method. For that use the 'network.method=\"NINA\"' option.\n\n3.- Represent your original distance matrix using the 'zero' method. For that use the 'network.method=\"zero\"' option.")

## GETTING NETWORKS ###
G<-graph.adjacency(dis2)
A<-as.network.matrix(dis2)
#######
#
if(save.distance==TRUE) write.table(file=save.distance.name,dis)
#
### BEGIN PLOT

Links<-as.matrix.network(A)
    png(filename="kk.png")
	vertis<-plot.network(A)
	file.remove("kk.png")
    dev.off()

if(is.na(label[1])==F & is.na(label.sub.str[1])==F)
{print("Multiple definition of labels")
label<-rep("",nrow(dis))}

if(is.na(label[1]) & is.na(label.sub.str[1])==F)
label<-substr(colnames(dis),label.sub.str[1],label.sub.str[2])

if(is.na(label[1])==F & is.na(label.sub.str[1]))
label<-label

if(is.na(label[1]) & is.na(label.sub.str[1]))
label<-colnames(dis)

Xmin<-min(vertis[,1])-(max(vertis[,1])-min(vertis[,1]))*0.1
Ymin<-min(vertis[,2])-(max(vertis[,2])-min(vertis[,2]))*0.1
Xmax<-max(vertis[,1])+(max(vertis[,1])-min(vertis[,1]))*0.1
Ymax<-max(vertis[,2])+(max(vertis[,2])-min(vertis[,2]))*0.1


## Colors ###

#if(modules==FALSE & is.na(bgcol[1]))
##bgcol<-colors()[sample(c(1,23,25:152,203:259,361:657),ncol(dis))]
#bgcol<-colour.scheme(def=bgcol,N=ncol(dis2))

if(modules==T)
	{
	comuni<-walktrap.community(G)
	tab1<-matrix(nrow=nrow(dis2),ncol=2)
	tab1<-as.data.frame(tab1)
	tab1[,1]<-label
	tab1[,2]<-comuni$membership
	#colo<-colors()[sample(c(1,23,25:152,203:259,361:657),length(unique(tab1[,2])))]
	#if(is.character(moduleCol[1])==T)
	#colo<-moduleCol
	colo<-colour.scheme(def=moduleCol,N=length(unique(tab1[,2])))
	tab1[which(tab1[,2]==1),3]<-colo[1]
	if(length(unique(tab1[,2]))>1)
	for(i in 2:length(unique(tab1[,2])))
	tab1[which(tab1[,2]==i),3]<-colo[i]
	colnames(tab1)<-c("Node_label","Module","Node_colour")
	bgcol<-tab1[,3]
#	out[[3]]<-tab1
#	names(out)<-c("Summary","Estimated Percolation Threshold","Module")
	write.table(file=modFileName,tab1,quote=F,row.names=FALSE)
	print(tab1)
	}


###

plot(c(1,1),xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax),type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
if(main=="summary")
mtext(paste("Network method: ",network.method,"      Indel method: ",indel.method,"\nSubstitution model: ", substitution.model,"      Pairwise deletion= ",pairwise.deletion,"\nAlpha= ",alpha,"      Combination method: ",combination.method,sep=""),font=2)

if(length(bgcol==1)) bgcol<-rep(bgcol,ncol(Links))

#Links<-as.matrix.network(Links)
#for (L1 in 1:(ncol(Links)-1))
#for (L2 in (L1+1):ncol(Links))
#if (Links[L1,L2]==1)
unos<-which(Links==1,arr.ind=T)
unos<-unos[which(unos[,1]<unos[,2]),]
EW1<-c()
EW2<-c()
OUTindels<-matrix(nrow=nrow(InDist),ncol=ncol(InDist))
OUTsubsts<-matrix(nrow=nrow(InDist),ncol=ncol(InDist))

for(Iunos in 1:nrow(unos))
	{
	L1<-unos[Iunos,1]
	L2<-unos[Iunos,2]

OUTindels[L1,L2]<-InDist[L1,L2]

# 3-take linked nodes and divide edge in the number of mutations it has

SuDist2<-as.matrix(dist.dna(x=alignUnique,model="N",pairwise.deletion=TRUE))

	InDisCor<-ceiling(InDist[L1,L2]/InScale)
	SuDisCor<-ceiling(SuDist2[L1,L2]/SuScale)

OUTsubsts[L1,L2]<-SuDist2[L1,L2]

LastEdgeIn<-InDist[L1,L2]-floor(InDist[L1,L2]/InScale)*InScale
LastEdgeWidthIn<-LastEdgeIn/InScale
if(LastEdgeWidthIn!=0)
EdgeWidthIn<-c(rep(lwd.mut,floor(InDist[L1,L2]/InScale)),LastEdgeWidthIn*lwd.mut)
else(EdgeWidthIn<-rep(lwd.mut,floor(InDist[L1,L2]/InScale)))

LastEdgeSu<-SuDist2[L1,L2]-floor(SuDist2[L1,L2]/SuScale)*SuScale
LastEdgeWidthSu<-LastEdgeSu/SuScale
if(LastEdgeWidthSu!=0)
EdgeWidthSu<-c(rep(lwd.mut,floor(SuDist2[L1,L2]/SuScale)),LastEdgeWidthSu*lwd.mut)
else(EdgeWidthSu<-rep(lwd.mut,floor(SuDist2[L1,L2]/SuScale)))

if(SuDist2[L1,L2]==0)
{EdgeWidthSu<-0}

EdgeWidth<-c(EdgeWidthIn,EdgeWidthSu)

EW1<-c(EW1,EdgeWidthIn)
EW2<-c(EW2,EdgeWidthSu)

	MUTS<-InDisCor+SuDisCor+1

	colores<-c(rep(colInd,InDisCor),rep(colSust,SuDisCor))

	nodo1<-vertis[L1,]
	nodo2<-vertis[L2,]

	edgeX_step<-(nodo1[1]-nodo2[1])/MUTS
	edgeY_step<-(nodo1[2]-nodo2[2])/MUTS

	segments(x0=nodo1[1],x1=nodo2[1],y0=nodo1[2],y1=nodo2[2],lwd=lwd.edge)

	edge.slope<-(nodo1[2]-nodo2[2])/(nodo1[1]-nodo2[1])
	edge.intercept<-nodo2[2]-edge.slope * nodo2[1]

	segment.slope<--1/edge.slope
	segment.angle<-atan(segment.slope)#/(pi/180)
	names(segment.angle)<-""
	segment.length<-(max(Links[,1])-min(Links[,1]))*0.01*cex.mut

		for (L3 in 1:(MUTS-1))
		{
		segment.intercept<-(nodo2[2]+(L3*edgeY_step))-segment.slope*(nodo2[1]+(L3*edgeX_step))
		x0=(nodo2[1]+(L3*edgeX_step))-segment.length*cos(segment.angle)
		x1=(nodo2[1]+(L3*edgeX_step))+segment.length*cos(segment.angle)
		y0=(nodo2[2]+(L3*edgeY_step))-segment.length*sin(segment.angle)
		y1=(nodo2[2]+(L3*edgeY_step))+segment.length*sin(segment.angle)
#		lines(x=c(x0,x1),y=c(y0,y1),lwd=3*lwd.mut,col=colores[L3])
		lines(x=c(x0,x1),y=c(y0,y1),lwd=3*EdgeWidth[L3],col=colores[L3])
		}
	}
#### PLOT CIRCLES #######
	points(x=vertis[,1],y=vertis[,2],pch=21,cex=4*cex.vertex,bg=bgcol)
	text(x=vertis[,1],y=vertis[,2],label,cex=0.8*cex.label)
#### END PLOT CIRCLES ###

	if(is.na(legend))
		{
		if(InScale!=1 | SuScale!=1) legend<-TRUE
		else(legend<-FALSE)
		}
	if(legend==TRUE)
	legend(legend.pos,c(paste(InScale,"indels"),paste(SuScale,"substitutions")), lwd=c(3*max(EW1), 3*max(EW2)), seg.len=0.1,col=c(colInd, colSust),bty=legend.bty,inset=0)

	row.names(OUTindels)<-label
	row.names(OUTsubsts)<-label
	colnames(OUTindels)<-label
	colnames(OUTsubsts)<-label

OUT<-list()
OUT[[1]]<-OUTsubsts
OUT[[2]]<-OUTindels
names(OUT)<-c("Substitutions","Indels")

if(modules==TRUE)
	{
	OUT[[3]]<-tab1
	names(OUT)[3]<-"modules"
	}

	if(silent==FALSE)
	{
	print(OUT)
	}
OUT

}
