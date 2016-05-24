compareGaussNetworks <-
function (Network1,Network2,Names=c("Network1","Network2"),Colors=c("yellow","blue","green"),interactiveMode=T,RhoThreshold=0.4,PValueThreshold=0.05)
{
	AcceptedClasses <- c("SIMoNeNet","WGCNANet")
	if (!(class(Network1) %in% AcceptedClasses) | !(class(Network2) %in% AcceptedClasses)) {stop("Bad classes used")}
	if (length(Names)!=2 & !is.character(Names)) {stop("Bad names used")}
	MergedEdges<-mergeGaussEdges(Network1$Edges,Network2$Edges,Names)
	AllNames<-c(Names,"Common")
	
	interactivePrompt<-"Press return for next plot..."
	
	Pos1<-which(MergedEdges$Network==Names[1])
	Pos2<-which(MergedEdges$Network==Names[2])
	CommonPos<-which(MergedEdges$Network=="Common")
	
	Area1<-length(which(MergedEdges$Network %in% c(Names[1],"Common")))
	Area2<-length(which(MergedEdges$Network %in% c(Names[2],"Common")))
	Middle<-length(which(MergedEdges$Network=="Common"))
	if (requireNamespace("VennDiagram",quietly=TRUE)) {Venn=VennDiagram::draw.pairwise.venn(Area1,Area2,Middle,fill=Colors[1:2],category=Names,cex=2,cat.cex=2,mar=0.05)} else {stop("VennDiagram package must be installed to use this function")}
	
	Rhos<-abs(MergedEdges$Rho)
	PValues<--log(MergedEdges$P.Value,10)
	maxpval<-99
	finitepos<-which(is.finite(PValues))
	if(length(finitepos)>0){maxpval<-max(PValues[finitepos])}
	infinitepos<-which(!is.finite(PValues))
	if(length(infinitepos)>0){PValues[infinitepos]<-maxpval}
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	Rhos1<-Rhos[Pos1]
	Rhos2<-Rhos[Pos2]
	CommonRhos<-Rhos[CommonPos]
	AllRhos<-list(Rhos1,Rhos2,CommonRhos)
	names(AllRhos)<-AllNames
	boxplot(AllRhos,col=Colors,main="Absolute values of rho in different networks",ylab="abs(Rho)")
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	PValues1<-PValues[Pos1]
	PValues2<-PValues[Pos2]
	CommonPValues<-PValues[CommonPos]
	AllPValues<-list(PValues1,PValues2,CommonPValues)
	names(AllPValues)<-AllNames
	boxplot(AllPValues,col=Colors,main="Spearman p-values in different networks",ylab="-log10(P-values)")
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	Connec1<-mean(table(c(Network1$Edges$node1,Network1$Edges$node2)))
	Connec2<-mean(table(c(Network2$Edges$node1,Network2$Edges$node2)))
	CommonConnec<-mean(table(c(MergedEdges$node1[CommonPos],MergedEdges$node2[CommonPos])))
	Connecs<-c(Connec1,Connec2,CommonConnec)
	names(Connecs)<-AllNames
	mp=barplot(Connecs,beside=T,col=Colors,main="Nodes connectivities",ylab="mean(connectivities)")
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	ColorVectors<-rep("black",nrow(MergedEdges))
	ColorVectors[Pos1]<-Colors[1]
	ColorVectors[Pos2]<-Colors[2]
	ColorVectors[CommonPos]<-Colors[3]
	plot(PValues~Rhos,pch=21,col="black",bg=ColorVectors,main="Spearman P-Values ~ Rhos",xlab="abs(Rho)",ylab="-log10(P-Value)")
	legend("topleft",legend=AllNames,pch=21,col="black",pt.bg=Colors,bg="white")
	abline(h=-log(PValueThreshold,10),lty=2,lwd=2)
	abline(v=RhoThreshold,lty=2,lwd=2)
	
	if(interactiveMode){line <- readline("Press return to close the plot...");dev.off()}
}
