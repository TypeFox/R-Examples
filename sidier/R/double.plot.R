double.plot <-
function(align=NA,indel.method="MCIC",substitution.model="raw", pairwise.deletion=TRUE,network.method.mut="percolation",network.method.pie="percolation",range=seq(0,1,0.01), addExtremes=FALSE,alpha.mut="info",alpha.pie="info",combination.method.mut="Corrected", combination.method.pie="Corrected",na.rm.row.col.mut=FALSE,na.rm.row.col.pie=FALSE, save.distance.mut=FALSE, save.distance.name.mut="DistanceMatrix_threshold_Mutations.txt", save.distance.pie=FALSE, save.distance.name.pie="DistanceMatrix_threshold_Pies.txt", modules=FALSE, modules.col=NA, bgcol=NA, label.col.mut="black", label.col.pie="black",label.mut=NA,label.pie=NA,label.sub.str.mut=NA,label.sub.str.pie=NA,colInd="red", colSust="black",lwd.mut=1,InScale=1, SuScale=1, lwd.edge=1.5,cex.mut=1,cex.label.mut=1,cex.label.pie=1,cex.vertex=1, main=c("Haplotypes","Populations"),NameIniPopulations=NA, NameEndPopulations=NA, NameIniHaplotypes=NA,NameEndHaplotypes=NA, cex.pie=1, HaplosNames=NA,offset.label=1.5)
{
if (modules==FALSE)
{
	if(is.na(bgcol[1]))
		{
		Aux<-GetHaplo(align=align, saveFile =FALSE, format = "fasta", seqsNames = NA,silent=T)
		if(is.na(bgcol))
		bgcol<-colour.scheme(def=bgcol,N=length(dimnames(Aux)[[1]]))
		col.pie<-NA
		}

	if(is.na(bgcol[1])==FALSE)
		bgcol->col.pie
}

if(length(main)==1 & main[1]=="summary")
{
ifelse(network.method.mut==network.method.pie,
NMH<-network.method.mut,NMH<-paste(network.method.mut,"&",network.method.pie))

ifelse(combination.method.mut==combination.method.pie,
CMH<-combination.method.mut,CMH<-paste(combination.method.mut,"&",combination.method.pie))

ifelse(alpha.mut==alpha.pie,
AMH<-alpha.mut,AMH<-paste(alpha.mut,"&",alpha.pie))

main[1]<-paste("Network method: ",NMH,";      Indel method: ",indel.method,";\nSubstitution model: ", substitution.model,";      Pairwise deletion= ",pairwise.deletion,";\nAlpha= ",AMH,";      Combination method: ",CMH,sep="")
main[2]<-paste("Haplotypes network (left) \nPopulations network (right)\n")
}

layout(matrix(c(1:2),ncol=2))

a<-mutation.network (align=align,indel.method=indel.method,substitution.model=substitution.model,pairwise.deletion=pairwise.deletion, network.method=network.method.mut,range=range,addExtremes=addExtremes,alpha=alpha.mut,combination.method=combination.method.mut, bgcol=bgcol, label.col=label.col.mut,modules=modules, moduleCol= modules.col,colInd=colInd, colSust=colSust,lwd.mut=lwd.mut,lwd.edge=lwd.edge,cex.mut=cex.mut,label=label.mut, cex.label=cex.label.mut,cex.vertex =cex.vertex, na.rm.row.col=na.rm.row.col.mut,save.distance=save.distance.mut, save.distance.name=save.distance.name.mut ,label.sub.str=label.sub.str.mut,silent=TRUE,InScale=InScale, SuScale=SuScale)
mtext(paste(main[1],"\n"),font=2)

if(modules==TRUE)
#col.pie<-unique(a[[3]][,3])
col.pie<-a[[3]][,3]

pie.network (align=align,indel.method=indel.method,substitution.model=substitution.model,pairwise.deletion=pairwise.deletion, network.method=network.method.pie,range=range,addExtremes=addExtremes,alpha=alpha.pie,combination.method=combination.method.pie, label=label.pie,label.col=label.col.mut,cex.label=cex.label.pie ,label.sub.str=label.sub.str.pie, na.rm.row.col=na.rm.row.col.pie, NameIniPopulations=NameIniPopulations, NameEndPopulations=NameEndPopulations, NameIniHaplotypes=NameIniHaplotypes, NameEndHaplotypes=NameEndHaplotypes,save.distance=save.distance.pie,save.distance.name=save.distance.name.pie, col.pie=col.pie,cex.pie=cex.pie,HaplosNames=HaplosNames,offset.label=offset.label)
mtext(paste(main[2],"\n"),font=2)
}
