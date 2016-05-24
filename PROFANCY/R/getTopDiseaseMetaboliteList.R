
##initialize data
initializeData<-function(){
   utils::data("envData",package="PROFANCY")
}
##initializeData
if(!exists("envData")) initializeData()


#### RandomWalk on graph
RandomWalk2igraph<-function(igraphM,VertexWeight,EdgeWeight=TRUE,gamma=0.7){
if(EdgeWeight==TRUE)
{
 adjM<-get.adjacency(igraphM,attr="weight") # convert igraph object to a weight matrix
}
if(EdgeWeight==FALSE)
{
 adjM<-get.adjacency(igraphM) # convert igraph object to a conventional matrix
}
adjM<-as.matrix(adjM)  #new
res<-rw(adjM,VertexWeight,gamma)
return(drop(res))
}

rw<-function(W,p0,gamma) {
   
   p0<-t(p0)
   p0 <- p0/sum(p0)
   PT <- p0
   
   k <- 0
   delta <- 1
   W<-as.matrix(W) #new
  Ng <- dim(W)[2]
  for (i in 1:Ng) {
      sumr<-sum(W[i,])
      if(sumr==0)
      {
      W[i,] <-numeric(length=length(W[i,]))
      }
      if(sumr>0)
      {
      W[i,] <- W[i,]/sum(W[i,])
      }
    }
   W <- t(W)
   
   while(delta>1e-10) {
      PT1 <- (1-gamma)*W
      PT2 <- PT1 %*% t(PT)
      PT3 <- (gamma*p0)
      PT4 <- t(PT2) + PT3
      delta <- sum(abs(PT4 - PT))
      PT <- PT4
      k <- k + 1
    }
    PT<-t(PT)
    return(PT)
}








##get network
##如果是用加通路信息的方法先得到加通路信息的网络 。如果不是，直接获得相应网络
##用户可以选择提供的两种网络，也可以自定义网络。如果自定义网络，输入的信息为两列，每一行是一个互作对。
getGraph<-function(network){
	if(network=="KEGG"){Net<-get("KEGGAddPathInfNetwork",envir=envData)}
	if(network=="EHMN"){Net<-get("EHMNAddPathInfNetwork",envir=envData)}
	return(Net)
}

##get seed nodes 
#如果疾病选择的是我们的疾病信息里存在的，那么根据他选择网络不通，分别得到不通seed。
#如果疾病是自定义用户输入的，那么seed也由用户自定义输入（KEGG Id），将其输入与不同网络中的代谢子取交集
getSeed<-function(diseaseName,network,seed,seedDefault){ 
if(seedDefault==TRUE){
	DiseaseInfList<-get("DiseaseInfList",envir=envData)
	ourDisease<-sapply(DiseaseInfList,function(x){x$OMIMName})
	if(diseaseName%in%ourDisease){
		location<-which(ourDisease%in%diseaseName)
		if(network=="KEGG"){
			seeds<-DiseaseInfList[[location]]$KEGGNetSeed
		}
		if(network=="EHMN"){
			seeds<-DiseaseInfList[[location]]$EHMNNetSeed
		}
	}
}
if(seedDefault==FALSE){
	KEGGAddPathInfNetwork<-get("KEGGAddPathInfNetwork",envir=envData)
	EHMNAddPathInfNetwork<-get("EHMNAddPathInfNetwork",envir=envData)
	if(network=="KEGG"){seeds<-intersect(seed,V(KEGGAddPathInfNetwork)$name)}
	if(network=="EHMN"){seeds<-intersect(seed,V(EHMNAddPathInfNetwork)$name)}

}
return(seeds)
}
####候选列表
getCandidates<-function(diseaseName,candidates,network,candidateDefault,seed,seedDefault){
KEGGAddPathInfNetwork<-get("KEGGAddPathInfNetwork",envir=envData)
EHMNAddPathInfNetwork<-get("EHMNAddPathInfNetwork",envir=envData)
if(candidateDefault==TRUE){
	SeedNode<-getSeed(diseaseName,network,seed,seedDefault) 
	a<-V(KEGGAddPathInfNetwork)$name
	KEGGNode<-a[which(substr(a,1,1)!="p")]
	b<-V(EHMNAddPathInfNetwork)$name
	EHMNNode<-b[which(substr(b,1,1)!="p")]
	if(network=="KEGG"){candidates<-setdiff(KEGGNode,SeedNode)}
	if(network=="EHMN"){candidates<-setdiff(EHMNNode,SeedNode)}
}
if(candidateDefault==FALSE){
	if(network=="KEGG"){candidates<-intersect(candidates,V(KEGGAddPathInfNetwork)$name)}
	if(network=="EHMN"){candidates<-intersect(candidates,V(EHMNAddPathInfNetwork)$name)}
}
	return(candidates)
}


##getPrioritization list
getPriList<-function(Graph,seed){
	library(igraph)
	V(Graph)$name->Vnames
	Attrlist<-list()
	VertexWeight<-rep(0,length=vcount(Graph))
	VertexWeight[which(Vnames%in%seed)]<-1
	visProbs <- RandomWalk2igraph(Graph,VertexWeight=VertexWeight,gamma=0.7,EdgeWeight=FALSE)
	names(visProbs)->Names
	which(substring(Names,1,4)%in%"path")->pathNodeIndex ##获得通路节点的索引
	#print(pathNodeIndex)
	if(length(pathNodeIndex)!=0){
		AllMetabolitesProbs<-visProbs[-pathNodeIndex] 			##new   所有代谢子的概率值
	}else{AllMetabolitesProbs<-visProbs}
	r1<-which(names(AllMetabolitesProbs) %in% seed)
	visProbs2<-AllMetabolitesProbs[-r1]
	AllMProbsNoSeedNode<-visProbs2[order(visProbs2,decreasing=TRUE)] ##new   所有代谢子不包括种子节点的概率值
	Attrlist[[1]]<-AllMetabolitesProbs
	Attrlist[[2]]<-AllMProbsNoSeedNode
	names(Attrlist)<-c("AllMetabolitesProbs","AllMProbsNoSeedNode")
	return(Attrlist)

}

#getTopDiseaseMetabolites<-function(diseaseName,network=c("KEGG","EHMN"),seed,candidates,
#        seedDefault=TRUE,candidateDefault=TRUE,showTop=30){

getTopDiseaseMetabolites<-function(diseaseName=NULL,network=c("KEGG"),seed=NULL,
	candidates=NULL,seedDefault=TRUE,candidateDefault=TRUE,showTop=30){

		if(is.null(diseaseName)==TRUE)diseaseName<-NA
		Graph<-getGraph(network)   ##get network
		SeedNode<-getSeed(diseaseName,network,seed,seedDefault)  ##get seed
		if(length(SeedNode)==0){print("The seeds you input are not in our network. Error!Pplease re-enter.");break}
		DiseaseNameInf<-paste("Disease Name is:" ,diseaseName)
		if(is.null(seed)==TRUE){
			SeedNodeInf1<-paste("The number of the seeds you input are: 0")			
		}else{SeedNodeInf1<-paste("The number of the seeds you input are: " ,length(seed))}
		SeedNodeInf2<-paste("The number of the seeds used in prioritizing the candidate metabolites are: " ,length(SeedNode))
		SeedNodeInf3<-paste("The seeds used in prioritizing the candidate metabolites are: ",paste(SeedNode,collapse=";"))
		NetworkInf<-paste("The network you used is :", network)
		Candidates<-getCandidates(diseaseName,candidates,network,candidateDefault,seed,seedDefault)
		if(length(Candidates)==0){print("The candidates you input are not in our network. Error!Pplease re-enter.");break}
		if(is.null(candidates)==FALSE){
			CandidatesInf1<-paste("The number of the candidate metabolites you input are: " ,length(candidates))
		}else{CandidatesInf1<-paste("The number of the candidate metabolites you input are: " ,0)}
		CandidatesInf2<-paste("The number of the candidate metabolites are prioritizing in this method are: " ,length(Candidates))
		DiseaseInf<-rbind(DiseaseNameInf,NetworkInf,SeedNodeInf1,SeedNodeInf2,SeedNodeInf3,CandidatesInf1,CandidatesInf2)
		DiseaseInf<-unname(DiseaseInf)
		print(DiseaseInf)
		#write.table(DiseaseInf,"DiseaseInf.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
		
		Attrist<-getPriList(Graph,SeedNode)	
		topMetabolitesID<-c()
		topMetabolitesPvalue<-c()
		if(candidateDefault==TRUE){
			topMetabolitesID<-names(Attrist$AllMProbsNoSeedNode)
			topMetabolitesPvalue<-Attrist$AllMProbsNoSeedNode
			show<-showTop
		}
		if(candidateDefault==FALSE){
			loci<-which(names(Attrist$AllMetabolitesProbs)%in%Candidates)
			if(length(loci)<=showTop)show<-length(loci)
			if(length(loci)>showTop)show<-showTop
			candidateInf<-Attrist$AllMetabolitesProbs[loci]
			orderedCandidateInf<-candidateInf[order(candidateInf,decreasing=TRUE)]
			topMetabolitesID<-names(orderedCandidateInf)
			topMetabolitesPvalue<-orderedCandidateInf
			
		}
		Rank<-c(1:length(topMetabolitesID))
		MetaboliteInf<-get("MetaboliteInf",envir=envData)
		topMetabolitesName<-sapply(topMetabolitesID,function(x){as.character(MetaboliteInf[which(MetaboliteInf[[1]]%in%x),"MetaboliteName"][[1]])})
		TopMetaboliteInf<-data.frame(Rank,topMetabolitesID,topMetabolitesName,topMetabolitesPvalue)
		names(TopMetaboliteInf)<-c("Rank","KEGGID","MetaboliteName","Score")
		TopMetaboliteInf<-TopMetaboliteInf[1:show,]		
		#write.table(TopMetaboliteInf,"TopMetaboliteInf.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
		return(TopMetaboliteInf)
}


##provided disease we studied.
getProvidedDiseaseName<-function(){
	DiseaseInfList<-get("DiseaseInfList",envir=envData)
	ProvidedDiseaseName<-rovidedDiseaseName<-sapply(DiseaseInfList,function(x){unlist(x$OMIMName)})
	return(ProvidedDiseaseName)

}

