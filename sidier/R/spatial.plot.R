spatial.plot <-
function(dis=NULL, align=NA, X=NULL, Y=NULL, indel.method="MCIC", substitution.model="raw",
pairwise.deletion=TRUE, alpha="info", combination.method="Corrected",
na.rm.row.col=FALSE, addExtremes=FALSE, NameIniPopulations=NA, NameEndPopulations=NA,
NameIniHaplotypes=NA, NameEndHaplotypes=NA, HaplosNames=NA, save.distance=FALSE,
save.distance.name="DistanceMatrix_threshold.txt", network.method="percolation",
range=seq(0,1,0.01), modules=FALSE, moduleCol=NA, modFileName="Modules_summary.txt",
bgcol="white", label.col="black", label=NA, label.sub.str=NA, label.pos= "b",
cex.label=1,cex.vertex=1,vertex.size="equal", plot.edges=TRUE, lwd.edge=1,to.ggmap=FALSE,
plot.ggmap=FALSE, zoom.ggmap=6, maptype.ggmap="satellite", label.size.ggmap=3)
{

	if(is.null(dis)==FALSE)
	pieSize<-rep(1,nrow(dis))

	if(is.null(dis))
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
		InDist<-MCIC(align=alignUnique, saveFile = F,silent=TRUE)
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
	#if(length(which(dis==0))!=nrow(dis) & is.na(zeros))
	#stop("\n\nSome of the off-diagonal elements in your matrix are zero. To display this matrix as a network, use the option 'zeros=\"min\"')")

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
	## ESTIMATING POPULATION DISTANCES FROM HAPLOTYPE DISTANCES AND HAPLOTYPE COMPOSITION:

		HaplosAll<-FindHaplo(align=align,saveFile=F) # That give new names to haplotypes
	
		if(is.na(NameIniHaplotypes)==F & is.na(NameEndHaplotypes)==F & is.na(HaplosNames))
		HaplosAll[,2]<-substr(HaplosAll[,1],NameIniHaplotypes,NameEndHaplotypes) #That will maintain the original names of haplotypes
		HaplosPop<-HapPerPop(saveFile=T,input=HaplosAll,NameIniPopulations=NameIniPopulations, NameEndPopulations=NameEndPopulations)

		if(is.na(HaplosNames)==F)
			{
			names.ori<-unique(HaplosAll[,2])
			names.fin<-matrix(nrow=nrow(HaplosAll))
			if(length(names.ori)!=length(HaplosNames)) stop("Incorrect number of haplotype names!")
			for (n1 in 1:length(names.ori))
				names.fin[which(HaplosAll[,2]==names.ori[n1]),]<-HaplosNames[n1]
			HaplosAll[,2]<-names.fin
			HaplosPop<-HapPerPop(saveFile=T,input=HaplosAll,NameIniPopulations=NameIniPopulations, NameEndPopulations=NameEndPopulations)
			}
	
		if(is.na(NameIniHaplotypes) & is.na(NameEndHaplotypes))
			{
			NameIniPopulations<-1
			if(length(which(as.matrix((strsplit(HaplosAll[,1],"")[[1]]))=="_"))==0) stop("Error in population names. It is recommended to use equal length sequence names with population and individual names separated by '_' (e.g., Pop01_id001...Pop23_id107). See ?pie.network for details.")
			if(is.na(NameEndPopulations)) warning("Population names defined by algorithm between character 1 and the first symbol '_' in sequences name.")
			NameEndPopulations<-which(as.matrix((strsplit(HaplosAll[,1],"")[[1]]))=="_")-1
			NameIniHaplotypes<-NameEndPopulations+1
			HaplosAll[,1]<-paste(HaplosAll[,1],HaplosAll[,2],sep="_")
			NameEndHaplotypes<-nchar(HaplosAll[1,1])
			HaplosAll[,1]<-HaplosAll[,1]
			colnames(dis)<-HaplosAll[match(colnames(dis),substr(HaplosAll[,1],1,nchar(colnames(dis)[1]))),2]		
			row.names(dis)<-colnames(dis)
			HaplosPop<-HapPerPop(saveFile=T,input=HaplosAll,NameIniPopulations=NameIniPopulations, NameEndPopulations=NameEndPopulations)
			NameIniHaplotypes<-1
			NameEndHaplotypes<-nchar(HaplosAll[1,2])
			}


		dis<-pop.dist(distances=dis,Haplos=HaplosPop[[1]], logfile=FALSE,saveFile=FALSE,NameIniHaplotypes=NameIniHaplotypes, NameEndHaplotypes=NameEndHaplotypes,NameIniPopulations=NameIniPopulations,NameEndPopulations=NameEndPopulations)

### vertex.size

	HP<-as.matrix(HaplosPop[[1]])
	pieSize<-rowSums(HP)/max(rowSums(HP)) # radius
	if(vertex.size=="area")
	pieSize<-sqrt(pieSize/pi)/max(sqrt(pieSize/pi)) # area
	if(vertex.size=="equal")
	pieSize<-rep(1,nrow(HP))
	if(vertex.size!="radius" & vertex.size!="area" & vertex.size!="equal")
	stop("wrong pie.size defined")

	}
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
#		print(paste("Threshold value:",j,"  Range to test: from ",min(range)," to ",max(range),sep=""))

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

#		if(is.null(colnames(dis)))
#		label<-c(1:ncol(dis))
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
#		print(paste("Threshold value:",j,"  Range to test: from ",min(range)," to ",max(range),sep=""))

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
		}
	j<-0
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
	j<-0
	}


### END ZERO THRESHOLD ###
#
#
### END THRESHOLD ESTIMATION ###
#
#
## WARNING IF percolation threshold is not found:
	if(is.na(j) & length(which(dis==0))!=nrow(dis))
	warning("\n\nPercolation threshold can not be estimated and some of the off-diagonal elements in your matrix are zero. Your distance matrix seems to provide low resolution. You may:\n\n1.- Redefine populations by meging those showing distance values of 0 before percolation threshold estimation. For that use the 'merge=TRUE' option \n\n2.- Represent your original distance matrix using the 'No Isolated Nodes Allowed' method. For that use the 'network.method=\"NINA\"' option.\n\n3.- Represent your original distance matrix using the 'zero' method. For that use the 'network.method=\"zero\"' option.")
#
#
#
## GETTING NETWORKS ###
G<-graph.adjacency(dis2)
A<-as.network.matrix(dis2)
#######
#
if(save.distance==TRUE) write.table(file=save.distance.name,dis)
#
### MODULES ESTIMATION

if(modules==TRUE)
	{
	comuni<-walktrap.community(G)
	tab1<-matrix(nrow=nrow(dis2),ncol=2)
	tab1<-as.data.frame(tab1)
	tab1[,1]<-row.names(dis2) #labels
	tab1[,2]<-comuni$membership
	colores<-colour.scheme(def=moduleCol,N=length(unique(tab1[,2])))

		comuni<-walktrap.community(G)
		tab1<-matrix(nrow=nrow(dis2),ncol=2)
		tab1<-as.data.frame(tab1)
		tab1[,1]<-label
		tab1[,2]<-comuni$membership
		colores<-tab1[,2]
#		bgcol<-colores
		colo<-colour.scheme(def=moduleCol,N=length(unique(tab1[,2])))
		if(is.character(moduleCol[1])==T)
		colo<-moduleCol
		tab1[which(tab1[,2]==1),3]<-colo[1]
		if(length(unique(tab1[,2]))>1)
		for(i in 2:length(unique(tab1[,2])))
		tab1[which(tab1[,2]==i),3]<-colo[i]
		colnames(tab1)<-c("Node_label","Module","Node_colour")
		bgcol<-tab1[,3]
		colores<-bgcol
		write.table(file=modFileName,tab1,quote=F,row.names=FALSE)


	}
if(modules==FALSE)
	colores<-bgcol


#### label names

if(is.na(label[1])==F & is.na(label.sub.str[1])==F)
{print("Multiple definition of labels")
label<-rep("",nrow(dis))}

if(is.na(label[1]) & is.na(label.sub.str[1])==F)
label<-substr(colnames(dis),label.sub.str[1],label.sub.str[2])

if(is.na(label[1])==F & is.na(label.sub.str[1]))
label<-label

if(is.na(label[1]) & is.na(label.sub.str[1]))
label<-colnames(dis)

#### label position

if(label.pos== "c"|label.pos== "center")
POS<-NULL
if(label.pos== "a"|label.pos== "above")
POS<-3
if(label.pos== "b"|label.pos== "below")
POS<-1
if(label.pos== "l"|label.pos== "left")
POS<-2
if(label.pos== "r"|label.pos== "right")
POS<-4

### BEGIN PLOT

if (is.null(X)|is.null(Y))
{
warning("Geographic coordinates not provided. Populations plotted according to the Fruchterman-Reingold algorithm.")
vertis1<-plot.network(A)
X<-vertis1[,1]
Y<-vertis1[,2]
}

rangeX<-max(X)-min(X)
rangeY<-max(Y)-min(Y)

if(rangeX == rangeY)
		{
		Xmin<-min(X)
		Xmax<-max(X)
		Ymin<-min(Y)
		Ymax<-max(Y)
		}

	if(rangeX < rangeY)
		{
		Xmin<-min(X)-(rangeY-rangeX)/2
		Xmax<-max(X)+(rangeY-rangeX)/2
		Ymin<-min(Y)
		Ymax<-max(Y)
		}

	if(rangeX > rangeY)
		{
		Ymin<-min(Y)-(rangeX-rangeY)/2
		Ymax<-max(Y)+(rangeX-rangeY)/2
		Xmin<-min(X)
		Xmax<-max(X)
		}

plot(X,Y,pch=21,bg=colores,xlab="",ylab="",axes=FALSE,xlim=c(0.99*Xmin,Xmax*1.01),ylim=c(0.99*Ymin,Ymax*1.01),cex=cex.vertex,col="black")

	if(plot.edges==TRUE)
	{
	vertis<-cbind(X,Y)
	Links<-as.matrix(A)
	for (L1 in 1:(ncol(Links)-1))
	for (L2 in (L1+1):ncol(Links))
	if (Links[L1,L2]==1)
		{
		nodo1<-vertis[L1,]
		nodo2<-vertis[L2,]
		segments(x0=nodo1[1],x1=nodo2[1],y0=nodo1[2],y1=nodo2[2],lwd=lwd.edge)
		}
	}

points(X,Y,pch=21,bg=colores,cex=cex.vertex*pieSize,col="black")
text(X,Y,label=label,pos=POS,cex=cex.label)

unos<-which(dis2==1,arr.ind=T)
links<-unos[which(unos[,1]>unos[,2]),]
row.names(links)<-paste("link",1:nrow(links))

colnames(dis2)<-row.names(dis2)

	out<-list()
	length(out)<-5
	names(out)<-c("location","colours","coordinates","network","links")
	out$location<-c(mean(X),mean(Y))
	out$colours<-colores
	out$coordinates<-data.frame(X,Y)
	out$network<-dis2
	out$links<-links
	names(out$coordinates) <- c("lon", "lat")


if(plot.ggmap==TRUE)
	{
	dev.new()
	hdf <- get_map(location=out$location,zoom=zoom.ggmap,maptype =maptype.ggmap)
	points <- out$coordinates

	if(plot.edges==TRUE)
		plot( 
		ggmap(hdf)
		+ annotate("segment", x = out$coordinates[links[,1],1], xend = out$coordinates[links[,2],1], y = out$coordinates[links[,1],2], yend = out$coordinates[links[,2],2], colour = "black")
		+ annotate("point",x=out$coordinates[,1], y=out$coordinates[,2],  ymin = 1, ymax = 2, colour = bgcol, size = 3.2*pieSize*cex.vertex)
		+ annotate('text',x=out$coordinates[,1], y=out$coordinates[,2], label = label, colour = label.col, size = 2.5*pieSize*cex.vertex)
			)

	if(plot.edges==FALSE)
		plot(
		ggmap(hdf)
		+ annotate("point",x=out$coordinates[,1], y=out$coordinates[,2],  ymin = 1, ymax = 2, colour = bgcol, size = 3.2*pieSize*cex.vertex)
		+ annotate('text',x=out$coordinates[,1], y=out$coordinates[,2], label = label, colour = label.col, size = 2.5*pieSize*cex.vertex)
			)
	}

colnames(dis2)<-row.names(dis2)
if(to.ggmap==TRUE)
	{
	print(out)
	out
	}

}
