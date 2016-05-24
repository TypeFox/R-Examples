
pie.network <-
function(align=NA,indel.method="MCIC",substitution.model="raw",pairwise.deletion=TRUE,network.method="percolation",range=seq(0,1,0.01), addExtremes=FALSE,alpha="info",combination.method="Corrected",na.rm.row.col=FALSE,NameIniPopulations=NA, NameEndPopulations=NA, NameIniHaplotypes=NA,NameEndHaplotypes=NA,save.distance=FALSE, save.distance.name="DistanceMatrix_threshold.txt", col.pie=NA, label.col="black",label=NA,label.sub.str=NA,cex.label=1,cex.pie=1,main="", HaplosNames=NA,offset.label=1.5,pie.size="equal")
{
 #require (igraph)
 #require (network)
#require(gridBase)
#require(grid)

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
	{	
		if(substr(HaplosAll[,1],NameIniHaplotypes,NameEndHaplotypes)[1]!="")
		HaplosAll[,2]<-substr(HaplosAll[,1],NameIniHaplotypes,NameEndHaplotypes) #That will maintain the original names of haplotypes
		if(substr(HaplosAll[,1],NameIniHaplotypes,NameEndHaplotypes)[1]=="")
		stop(paste("Wrong haplotype names!  According to your input, haplotype names must be contained in sequence names between position",NameIniHaplotypes,"and",NameEndHaplotypes,", but it is not the case!"))
	HaplosPop<-HapPerPop(saveFile=T,input=HaplosAll,NameIniPopulations=NameIniPopulations, NameEndPopulations=NameEndPopulations)
	}

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
		if(is.na(NameIniPopulations)==FALSE & is.na(NameEndPopulations)==FALSE)
			{
#			if(length(which(as.matrix((strsplit(HaplosAll[,1],"")[[1]]))=="_"))==0) stop("Error in population names. It is recommended to use equal length sequence names with population and individual names separated by '_' (e.g., Pop01_id001...Pop23_id107). See ?pie.network for details.")
			NameIniHaplotypes<-nchar(HaplosAll[1,1])+1
			HaplosAll[,1]<-paste(HaplosAll[,1],HaplosAll[,2],sep="_")
			NameEndHaplotypes<-nchar(HaplosAll[1,1])
			HaplosAll[,1]<-HaplosAll[,1]
			colnames(dis)<-HaplosAll[match(colnames(dis),substr(HaplosAll[,1],1,nchar(colnames(dis)[1]))),2]		
			row.names(dis)<-colnames(dis)
			HaplosPop<-HapPerPop(saveFile=T,input=HaplosAll,NameIniPopulations=NameIniPopulations, NameEndPopulations=NameEndPopulations)
			NameIniHaplotypes<-1
			NameEndHaplotypes<-nchar(HaplosAll[1,2])
			NameEndPopulations<-NameEndPopulations-NameIniPopulations+1
			NameIniPopulations<-1
			}
		if(is.na(NameIniPopulations) & is.na(NameEndPopulations))
			{
			NameIniPopulations<-1
			if(length(which(as.matrix((strsplit(HaplosAll[,1],"")[[1]]))=="_"))==0) stop("Error in population names. It is recommended to use equal length sequence names with population and individual names separated by '_' (e.g., Pop01_id001...Pop23_id107). See ?pie.network for details.")
			warning("Population names defined by algorithm between character 1 and the first symbol '_' in sequences name.")
			NameEndPopulations<-(which(as.matrix((strsplit(HaplosAll[,1],"")[[1]]))=="_")-1)[1]
#silen26-5	NameIniHaplotypes<-NameEndPopulations+1 #(cambio x abajo)
			NameIniHaplotypes<-(nchar(HaplosAll[1,1])+2)
			HaplosAll[,1]<-paste(HaplosAll[,1],HaplosAll[,2],sep="_")
			NameEndHaplotypes<-nchar(HaplosAll[1,1])
			HaplosAll[,1]<-HaplosAll[,1]
			colnames(dis)<-HaplosAll[match(colnames(dis),substr(HaplosAll[,1],1,nchar(colnames(dis)[1]))),2]		
			row.names(dis)<-colnames(dis)
			HaplosPop<-HapPerPop(saveFile=T,input=HaplosAll,NameIniPopulations=NameIniPopulations, NameEndPopulations=NameEndPopulations)
			NameIniHaplotypes<-1
			NameEndHaplotypes<-nchar(HaplosAll[1,2])
			}

		}


####### PROBAR A REDEFINIR LOS NOMBRES DE LAS SECUENCIAS NADA MÁS EMPEZAR Y VER SI ASÍ LO HACE BIEN...


	dis<-pop.dist(distances=dis,Haplos=HaplosPop[[1]], logfile=FALSE,saveFile=FALSE,NameIniHaplotypes=NameIniHaplotypes, NameEndHaplotypes=NameEndHaplotypes,NameIniPopulations=NameIniPopulations,NameEndPopulations=NameEndPopulations)
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
## GETTING NETWORKS ###
G<-graph.adjacency(dis2)
A<-as.network.matrix(dis2)
#######
#
if(save.distance==TRUE) write.table(file=save.distance.name,dis)
#
### BEGIN PLOT

if(is.na(col.pie[1])==FALSE)
if(length(col.pie)!=length(unique(HaplosAll[,2])))
col.pie<-NA

#col.pie<-colors()[sample(c(1,23,25:152,203:259,361:657),length(unique(HaplosAll[,2])))]
if(is.na(col.pie[1]))
col.pie<-colour.scheme(def=col.pie,N=length(unique(HaplosAll[,2])))

Links<-as.matrix.network(A)
#vertis<-plot.network(A)

if(is.na(label[1])==F & is.na(label.sub.str[1])==F)
{print("Multiple definition of labels")
label<-rep("",nrow(dis))}

if(is.na(label[1]) & is.na(label.sub.str[1])==F)
label<-substr(colnames(dis),label.sub.str[1],label.sub.str[2])

if(is.na(label[1])==F & is.na(label.sub.str[1]))
label<-label

if(is.na(label[1]) & is.na(label.sub.str[1]))
label<-colnames(dis)

### PLOT PIES ###
#    png(file="mygraphic.png",width=400,height=350)
#    plot(x=rnorm(10),y=rnorm(10),main="example")
#    dev.off()

#vertis<-plot.network(A,vertex.col=NULL,label=NULL,usearrows=0,vertex.cex=0.4,interactive=F, label.pos=5,label.col=NULL,label.cex=0.8*cex.label,edge.col="white",vertex.border="white")
vertis<-plot.network(A,vertex.col="white",label="",usearrows=0,vertex.cex=0.4,interactive=F, label.pos=5,label.cex=0.8*cex.label,edge.col="white",vertex.border="white")


if(main=="summary")
mtext(paste("Network method: ",network.method,"      Indel method: ",indel.method,"\nSubstitution model: ", substitution.model,"      Pairwise deletion= ",pairwise.deletion,"\nAlpha= ",alpha,"      Combination method: ",combination.method,sep=""),font=2)

#	library(lattice)
#	library(gridBase)
#	library(grid)
	
#	HaplosAll<-FindHaplo(align=align,saveFile=F) #silenciado 15mayo
	if(exists("HaplosNames"))
	HaplosAll[,2]<-HaplosNames
#	HaplosPop<-HapPerPop(saveFile=F,input=HaplosAll,NameIniPopulations=NameIniPopulations, NameEndPopulations=NameEndPopulations)#silenciado 15mayo
	print(as.data.frame(rbind(col.pie,HaplosPop[[1]])))

	oldpar <- par(no.readonly = TRUE)

	x<-vertis[,1]
	y<-vertis[,2]

	vps <- baseViewports()
	par(new = TRUE)
	pushViewport(vps$inner, vps$figure,vps$plot)
	
	maxpiesize <- unit(1, "inches")
	totals <- 1
	sizemult <- rep(0.5,nrow(dis2))

	HP<-as.matrix(HaplosPop[[1]])

	pushViewport(viewport(x = unit(x[1],"native"), y = unit(y[1],"native"), width = sizemult[1] *maxpiesize, height = sizemult[1] *maxpiesize))
	grid.rect(gp = gpar(col = "white",fill = NULL, lty = "blank"))
	par(plt = gridPLT(), new = TRUE)
	popViewport()
	popViewport(3)

	par(oldpar)

pieSize<-rowSums(HP)/max(rowSums(HP)) # radius
if(pie.size=="area")
pieSize<-sqrt(pieSize/pi)/max(sqrt(pieSize/pi)) # area
if(pie.size=="equal")
pieSize<-rep(1,nrow(HP))
if(pie.size!="radius" & pie.size!="area" & pie.size!="equal" & pie.size!="points")
stop("wrong pie.size defined")

cor<-(rowSums(HP)/max(rowSums(HP)))

	POS<-rep(3,nrow(vertis))
	for(i in 1:nrow(vertis))
#	if(vertis[i,2]==max(vertis[,2]))
#	POS[i]<-1
	text(x=vertis[,1],y=vertis[,2],label,cex=0.8*cex.label,pos=POS,offset=offset.label)
#	text(x=vertis[,1],y=(vertis[,2]-vertis[,2]*(1-pieSize)*0.009),label,cex=0.8*cex.label,pos=POS,offset=offset.label) Para que la distancia entre pie y main sea cte.


#Edges
#	pushViewport(viewport(x = unit(x[i],"native"), y = unit(y[i],"native"), width = sizemult[i] *maxpiesize, height = sizemult[i] *maxpiesize))
for (L1 in 1:(ncol(Links)-1))
for (L2 in (L1+1):ncol(Links))
if (Links[L1,L2]==1)
{
	nodo1<-vertis[L1,]
	nodo2<-vertis[L2,]

	edge.slope<-(nodo1[2]-nodo2[2])/(nodo1[1]-nodo2[1])
	edge.intercept<-nodo2[2]-edge.slope * nodo2[1]

	radio<-0.05
	Npart<-13

	edgeX_step1<-(nodo1[1]-nodo2[1])/Npart
	edgeY_step1<-(nodo1[2]-nodo2[2])/Npart
	edgeX_step2<-(Npart-1)*(nodo1[1]-nodo2[1])/Npart
	edgeY_step2<-(Npart-1)*(nodo1[2]-nodo2[2])/Npart


		x0=nodo1[1]
		y0=nodo1[2]
		x1=nodo2[1]
		y1=nodo2[2]

#		x0=(nodo2[1]+edgeX_step1)-segment.length*cos(segment.angle)
#		x1=(nodo2[1]+edgeX_step2)+segment.length*cos(segment.angle)
#		y0=(nodo2[2]+edgeY_step1)-segment.length*sin(segment.angle)
#		y1=(nodo2[2]+edgeY_step2)+segment.length*sin(segment.angle)
		lines(x=c(x0,x1),y=c(y0,y1))

#	segments(x0=x0,x1=x1,y0=y0,y1=y1,lwd=1*lwd.edge)
}

if(pie.size=="points")
	{
	if(nrow(HP)<37)
		{
		coloresAux<-c("blue","green2","red","yellow","DarkOrchid1","gray51","chocolate","cyan4","saddle brown","aquamarine","chartreuse","chocolate1","DarkOrchid3","gray18","gold","DarkOrchid4","green4","gray29", "sienna3","tan1","blue4","limegreen","gray73","bisque3","deeppink","red4","OliveDrab4","gray95", "salmon","DeepPink4","green yellow","gray4","hot pink","pink2","dark orange","gold3")
		col.pie<-coloresAux[1:nrow(HP)]
		}
	else
		{
		coloresAux<-colors()[sample(c(1,23,25:152,203:259,361:657),nrow(HP))]
		col.pie<-coloresAux[1:nrow(HP)]
		}
	points(vertis,cex=(1.5*cex.pie),bg=col.pie,pch=21)
	}


if(pie.size!="points")
	{

	oldpar <- par(no.readonly = TRUE)

	x<-vertis[,1]
	y<-vertis[,2]

	vps <- baseViewports()
	par(new = TRUE)
	pushViewport(vps$inner, vps$figure,vps$plot)
	
	maxpiesize <- unit(1, "inches")
	totals <- 1
	sizemult <- rep(0.5,nrow(dis2))

	HP<-as.matrix(HaplosPop[[1]])

	if(pie.size!="points")
		for (i in 1:nrow(HP)) 
		{
		pushViewport(viewport(x = unit(x[i],"native"), y = unit(y[i],"native"), width = sizemult[i] *maxpiesize, height = sizemult[i] *maxpiesize))
		grid.rect(gp = gpar(col = "white",fill = NULL, lty = "blank"))
		par(plt = gridPLT(), new = TRUE)
		pie(HP[i,],radius=(1*cex.pie*pieSize[i]),labels="",col=col.pie)
		popViewport()
		}
	popViewport(3)
	
	par(oldpar)
	}


}
