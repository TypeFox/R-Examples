#Required packages:
#library(GO.db)   
#library(cluster)    

#Preloaded functions:
rm.na<-function (x) 
{
    ind <- is.na(x) | is.nan(x) | is.infinite(x)
    return(x[!ind])
}

which.na<- function(x){ which(is.na(x)) }

match.list<-function(string, list){
	LenG=0
	NameG=character()
	for (i in 1:length(list)){
		LenG[i]=length(list[[i]])
		if (LenG[i]>0) {NameG=c(NameG, names(list[i]))}
	}
	CumG=cumsum(LenG)
	UcomG=unique(CumG)+1
	X=match(string, as.character(unlist(list)))
	Gind=NameG[findInterval(X, UcomG)+1]
	match(Gind, names(list))
}

lcombine<-function(N,m){
lfactorial(N)-lfactorial(N-m)-lfactorial(m)}

LP1<-function(N,n1,n2,m){
lcombine(N,m)+lcombine(N-m,n1-m)+lcombine(N-n1,n2-m)-lcombine(N,n1)-lcombine(N,n2)}



### Section I: calculate "Significant Protein Pairs", and make plot of protein clusters ###
#' @name SignificantPairs
#' @title Identify functionally associated protein pairs
#' @description This function uses the PAND distribution to calculate p-values (or probabilities) for each pair of proteins with at least one common neighbor in a protein-protein interaction network. It returns protein pairs with significant p-values (or probabilities).
#' @param PPIdb A two-column data frame consisting of binary interactions where each row represents an undirected edge (interaction) between two nodes (proteins) from two columns.  
#' @param Lambda Weight of direct interactions in the PAND algorithm. Lamda has different biological meanings with different values: "0" indicates that a direct link gives no information on the functional association; "1" indicates that a direct link is as informative as sharing one common neighbor (defined as an indirect link) on the functional association; "2" (or greater integer) indicates that a direct link is more informative than an indirect link. Since direct links should represent stronger evidence of functional associations than indirect links, we recommend using "2" as Lamda. 
#' @param pvalue logical; if TRUE, p-values for protein pairs will be calculated using PAND; if FALSE, probabilities will be calculated.
#' @usage SignificantPairs(PPIdb, Lambda=2, pvalue=FALSE)
#' @return This function returns a data frame with column names: "Sym_A", "Sym_B", "p_value" and "CommonNeighbor". "Sym_A" and "Sym_B" are a pair of nodes that share a significant functional linkage. "p_value" or "Probability" (calculated by the PAND algorithm) measures the significance of the linkage. "CommonNeighbor" is the number of shared nodes.
#' @export
#' @rdname SignificantPairs-methods
#### imported packages ###########
#' @import GO.db
#' @import cluster
#' @examples
#' ## not run 
#' ## data(dfPPI)
#' ## OrderAll=SignificantPairs(dfPPI)
#' @seealso \code{\link{ProteinCluster}}, \code{\link{KEGGpredict}}, \code{\link{GOpredict}}, \code{\link{SignificantSubcluster}}


SignificantPairs <- function(PPIdb, Lambda=2, pvalue=FALSE){
	DataBase=PPIdb
	proall=unique(rm.na(c(DataBase[,1],DataBase[,2])))
	u=length(proall)
	LenH=length(proall)
	proteinvector=vector(mode='list', length=u)
	Len=0
	for(i in 1:u){
		t1=DataBase[DataBase[,1]==proall[i],2]
		t2=DataBase[DataBase[,2]==proall[i],1] 
		t=unique(rm.na(c(t1,t2)))
		proteinvector[[i]]=t
		Len[i]=length(t)
	}
	names(proteinvector)<- proall
	Xs=mean(Len)
	XsS=mean(Len^2)
	
	ProA=vector(mode='character', length=6000000)
	ProB=vector(mode='character', length=6000000)
	Pr1=vector(mode='numeric', length=6000000)
	Pr3=vector(mode='numeric', length=6000000)
	Pr5=vector(mode='numeric', length=6000000)
	CmNeib=vector(mode='numeric', length=6000000)
	Direct=vector(mode='numeric', length=6000000)
	Id=0
  INDEX=0
  
  if (pvalue) {
    for (i in  1:(u-1)){
      n1=length(proteinvector[[i]])
      for (k in (i+1):u ){
        n2=length(proteinvector[[k]])
        inter=intersect(proteinvector[[i]],proteinvector[[k]])
        Dir=length(rm.na(match(proall[i], proteinvector[[k]])))*Lambda
        m=length(inter)
        INDEX=INDEX+1
        if(m!=0){
          Id=Id+1
          ProA[Id]=proall[i]
          ProB[Id]=proall[k]
          Pr1[Id]=log(sum(exp(LP1(u,n1,n2,m:min(n1,n2)))))   # for LP1, Phi is 1.
          Phi=sum(exp(LP1(u,n1,n2,0:min(n1,n2))+(0:min(n1,n2))*log(XsS)-(2*(0:min(n1,n2)))*log(Xs)))
          Pr3[Id]=log(sum(exp(LP1(u,n1,n2,m:min(n1,n2))+(m:min(n1,n2))*log(XsS)-(2*(m:min(n1,n2)))*log(Xs))))-log(Phi)
          Phi=sum(exp(LP1(u,n1+Dir,n2+Dir,0:min(n1+Dir,n2+Dir))+ (0:min(n1+Dir,n2+Dir))*log(XsS)-(2*(0:min(n1+Dir,n2+Dir)))*log(Xs)))
          Pr5[Id]=log(sum(exp(LP1(u,n1+Dir,n2+Dir,(m+Dir):min(n1+Dir,n2+Dir))+ ((m+Dir):min(n1+Dir,n2+Dir))*log(XsS)-(2*((m+Dir):min(n1+Dir,n2+Dir)))*log(Xs))))-log(Phi)
          CmNeib[Id]=m
          Direct[Id]=Dir/2
        }
        #if(is.element(INDEX/choose(u,2), 1:10/10)){
          #print(paste(100*INDEX/choose(u,2), '%', ' completed', sep=''))
        #} 
      }
    }
    DataHomoIVplus=data.frame(Sym_A=ProA, Sym_B=ProB, p_pvalue=Pr5, CommonNeighbor=CmNeib, stringsAsFactors=F)
  }
	else {
	  for (i in  1:(u-1)){
	    n1=length(proteinvector[[i]])
	    for (k in (i+1):u ){
	      n2=length(proteinvector[[k]])
	      inter=intersect(proteinvector[[i]],proteinvector[[k]])
	      Dir=length(rm.na(match(proall[i], proteinvector[[k]])))*Lambda
	      m=length(inter)
	      INDEX=INDEX+1
	      if(m!=0){
	        Id=Id+1
	        ProA[Id]=proall[i]
	        ProB[Id]=proall[k]
	        #Pr1[Id]=LP1(u,n1,n2,m)
	        #Pr3[Id]=LP1(u,n1,n2,m)+ log(XsS^m/Xs^(2*m))
	        Pr5[Id]=LP1(u,n1+Dir,n2+Dir,m+Dir)+ log(XsS^(m+Dir)/Xs^(2*(m+Dir)))
	        CmNeib[Id]=m
	        Direct[Id]=Dir/2
	      }
	      #if(is.element(INDEX/choose(u,2), 1:10/10)){
	        #print(paste(100*INDEX/choose(u,2), '%', ' completed', sep=''))	      
	      #}
	    }
	  }
	  DataHomoIVplus=data.frame(Sym_A=ProA, Sym_B=ProB, Probability=Pr5, CommonNeighbor=CmNeib, stringsAsFactors=F)
	}
	a=match('', DataHomoIVplus[,1])
	DataHomoIVplus=DataHomoIVplus[1:(a-1),]   #remove "0" lines
	Homo=DataHomoIVplus
	u=length(proall)
	CutOff=0.05/choose(u,2)
	HomoHM=Homo[order(Homo[,3]),]
	orderall=HomoHM[HomoHM[,3]< log(CutOff),]
	Pro=unique(c(orderall[,1],orderall[,2]))
	orderall
}
#OrderAll=SignificantPairs(dfPPI)


#' @title Cluster proteins based on significant protein pairs
#'
#' @description This function uses the p-values (or probabilities) derived from the PAND algorithm to perform agglomerative hierarchical clustering (using the unweighted group average) for proteins that form significant protein pairs. 
#'
#' @param Pfile A data frame returned from the function SignificantPairs()
#' @param Plot If FALSE, a dendrogram will NOT be generated
#' @param TextScaler Scale the size of the label in the generated PDF file
#' @param height The height of the generated PDF file
#' @param width The width of the generated PDF file
#' @usage ProteinCluster(Pfile, Plot=FALSE, TextScaler=50, height=10, width)
#' @return This function returns an object in the class "dendrogram". If the argument "Plot" is "TRUE", it will also plot the dendrogram.
#' @seealso \code{\link{SignificantPairs}}, \code{\link{KEGGpredict}}, \code{\link{GOpredict}}, \code{\link{SignificantSubcluster}}
#' @examples
#' ## not run
#' ## data(dfPPI)
#' ## OrderAll=SignificantPairs(dfPPI)
#' ## dendMap=ProteinCluster(Pfile=OrderAll, Plot=TRUE, TextScaler=30)
#' @export
#' @rdname ProteinCluster-methods

ProteinCluster<-function(Pfile, Plot=FALSE, TextScaler=50, height=10, width){ 
	orderup=Pfile
	names(orderup)[3] <- 'value'
	#attach(orderup)
	fn1=ecdf(orderup$value)
	score=numeric()
	a=0
	for (i in 1:dim(orderup)[1]){    
		a=a+1
		score[a]=fn1(orderup$value[i])
	}
	#detach(orderup)
	ProEx=unique(c(orderup[,1], orderup[,2]))

	ScoreData=data.frame(X=orderup[,1], Y=orderup[,2], score=score, stringsAsFactors=F) 
	rm(score)

	#attach(ScoreData)
	MatrixGraph=matrix(10,length(ProEx),length(ProEx),dimnames=list(ProEx, ProEx))
	X=ScoreData[,1]
  Y=ScoreData[,2]
  score=ScoreData[,3]
  for (i in 1:nrow(ScoreData)){
		MatrixGraph[X[i], Y[i]]= score[i]
		MatrixGraph[Y[i], X[i]]= score[i]
	}
	#detach(ScoreData)
	Map=agnes(MatrixGraph,diss=TRUE)

	hMap<- as.hclust(Map)
	dendMap<- as.dendrogram(hMap)
	nP <- list(col=3:2, cex=c(1.0, 0.5), pch= 21:22, bg= c("light blue", "pink"), lab.cex = 0.75, lab.col = "tomato")
 	addE <- function(n) {
		if(!is.leaf(n)) {
			attr(n, "edgePar") <- list(p.col="plum")
			attr(n, "edgetext") <- paste(attr(n,"members"),"members")
		}
		n
	}

	if (Plot==TRUE){
		if (missing(width)){L=min(length(ProEx)/10, 200)}
		#pdf(GraphFile, height, width=L)
		plot(dendMap, nodePar=list(lab.cex=L/TextScaler, cex=c(0,0), lab.col='tomato'))   # this is the second output.
		#dev.off()
	}
	dendMap
}
### dendMap=ProteinCluster(Pfile=OrderAll, Plot=TRUE, TextScaler=30)  # Generate a dendrogram of proteins.


### Section II: predict GO terms and KEGG terms (only) for proteins ###

#' @title Predict KEGG pathway annotations for proteins
#'
#' @description This function uses a direct annotation scheme to predict KEGG pathway annotations for proteins in the network derived with the PAND algorithm.
#'
#' @param Pfile A data frame returned from the function "SignificantPairs"
#' @param PPIdb A 2-column data frame consisting of binary interactions where each row i.e. c(A, B) represents an undirected edge (interaction) between gene A and gene B. 
#' @param Gene2Annotation A list that maps KEGG pathway ID to genes. The names should be gene symbols and the elements should be KEGG pathways. i.e. 
#'
#'	$IGSF5
#'
#'	[1] "hsa04530" "hsa05120"
#'
#'	$OR2T8
#'
#'	[1] "hsa04740"
#'
#'	$hCG_1776980
#'
#'	[1] "hsa00020" "hsa00190" "hsa01100" "hsa05010" "hsa05012" "hsa05016"
#'
#'	......
#'
#' @param p_value A cut-off for p-values from Fishers exact test when predicting KEGG pathway annotations
#' @param IDtoNAME A table that maps KEGG pathway ID to KEGG pathway names.
#' @usage KEGGpredict(Pfile, PPIdb, Gene2Annotation, p_value=0.001, IDtoNAME)
#' @return This function returns a data frame with column names: "Symbol", "KEGGID", "PathName", "Ratio" and "Pvalue". "Symbol" is the name of the 
#' node this function predicts KEGG annotation for. "KEGGID" and "PathName" are the predicted KEGG pathway identifier and the pathway name, respectively. 
#' "Ratio" is the proportion of neighboring nodes that have the predicted KEGG annotation. "Pvalue" is calculated from the "Ratio" by Fishers exact test. 
#' @seealso \code{\link{SignificantPairs}}, \code{\link{ProteinCluster}}, \code{\link{GOpredict}}, \code{\link{SignificantSubcluster}}
#' @export
#' @rdname KEGGpredict-methods
KEGGpredict<-function(Pfile, PPIdb, Gene2Annotation, p_value=0.001, IDtoNAME){
	#DataBase=read.delim(PPIdb, as.is=T)
	DataBase=PPIdb
	LenH=length(unique(rm.na(c(DataBase[,1],DataBase[,2]))))
	#keggvector=dget(Gene2Annotation)
	keggvector=Gene2Annotation
	orderup=Pfile
	proall=unique(c(orderup[,1], orderup[,2]))
	keggidall=list()
	neibor=list()
	aa=0
	for (i in proall){           
		aa=aa+1
		t1=orderup[orderup[,1]==i,2]
		t2=orderup[orderup[,2]==i,1] 
		t=unique(rm.na(c(t1,t2)))
		keggid=character()
		for (k in t ){
			id= unique(keggvector[[k]])
			keggid= c(keggid,id)
		}
		keggidall[[aa]]=keggid
		neibor[[aa]]=length(t)
	}
	names(keggidall) <- proall
	names(neibor)    <- proall
 
	allkegg=sort(table(unlist(keggvector)), decreasing=T)
	keggname=names(allkegg)
	keggcount=unname(allkegg)
	kgLen=length(keggvector)
	allKanno=sum(keggcount)
	KEGGannoAssign=character()
	ovlap=numeric()
	Pvalue=numeric()
	ind=0
	for (i in proall){
		ind=ind+1
		if (length(AA<-keggidall[[i]])>0){
			z=numeric()
			ov=numeric()
			for (k in unique(AA)){
				x1=c(keggcount[match(k, keggname)]-sum(AA==k), LenH-keggcount[match(k, keggname)]-(neibor[[i]]-sum(AA==k)))  # LenH is the number of human proteins in PPI data.
				x2=c(sum(AA==k), neibor[[i]]-sum(AA==k))
				xx=data.frame(x1=x1,x2=x2)
				z=c(z,fisher.test(xx)$p.value)
				ov=c(ov, sum(AA==k))
			}
			if (min(z)< p_value){
				KEGGannoAssign[ind]=unique(AA)[which.min(z)]
				Pvalue[ind]=min(z)
				ovlap[ind]=ov[which.min(z)]
			}
		}
	}

	TrueAnnoKEGG=rep(NA, length(proall))
	Ra=rep(NA, length(proall)) 
	Rao=rep(NA, length(proall)) 
	for (i in 1:length(KEGGannoAssign)){
		A<-KEGGannoAssign[i][!is.element(rm.na(KEGGannoAssign[i]),keggvector[[proall[i]]])]
		if (length(A)>0) {
			TrueAnnoKEGG[i]=A
			Ra[i]=paste(ovlap[i],'/', neibor[[i]], sep='')
			Rao[i]=ovlap[i]/neibor[[i]]
		}
	}

	Kp=IDtoNAME
	X=which(!is.na(TrueAnnoKEGG))
	keggPath=character()
	for (i in TrueAnnoKEGG[X]){
		pathID <- sub('hsa', '', i)
		pathName <- Kp[which(Kp[,1]==pathID),2]
		keggPath=c(keggPath, pathName)
	}
	#browser()
	KEGGpredict=data.frame(Symbol=proall[X], KEGGID=TrueAnnoKEGG[X], PathName=keggPath, Ratio=Ra[X], Pvalue=Pvalue[X])
	KEGGpredict
}
###KP=KEGGpredict(Pfile=OrderAll, PPIdb=dfPPI, Gene2Annotation='GENE2KEGG', p_value=0.001, IDtoNAME=KEGGID2NAME)


#' @title Predict GO annotations for proteins
#'
#' @description This function uses a direct annotation scheme to predict GO annotations for proteins in the network derived with the PAND algorithm
#'
#' @param Pfile A data frame returned from the function "SignificantPairs"
#' @param PPIdb A 2-column data frame consisting of binary interactions where each row i.e. c(A, B) represents an undirected edge (interaction) between gene A and gene B. 
#' @param Gene2Annotation A list that maps GO ID to genes. The names should be gene symbols and the elements should be GO IDs. i.e. 
#'
#'	$SHC1
#'
#'	[1] "GO:0005158" "GO:0005068" "GO:0005159" "GO:0070435"
#'
#'	$POU5F1
#'
#'	[1] "GO:0035413" "GO:0003130" "GO:0060391" "GO:0090308" "GO:0060965" "GO:0035198"
#'
#'	$FGF12
#'
#'	[1] "GO:0008201"
#'
#'	......
#'
#' @param p_value A cut-off for p-values from Fishers exact test when predicting GO annotations
#' @usage GOpredict(Pfile, PPIdb, Gene2Annotation, p_value=0.001)
#' @return This function returns a data frame with column names: "Symbol", "GOID", "GOterm", "Ratio" and "Pvalue". "Symbol" is the name of the node this function predicts GO annotation for. "GOID" and "GOterm" are the predicted GO ID and the GO term, respectively. "Ratio" is the proportion of neighboring nodes that have the predicted GO annotation. "Pvalue" is calculated from the "Ratio" by Fishers exact test. 
#' @seealso \code{\link{SignificantPairs}}, \code{\link{ProteinCluster}}, \code{\link{KEGGpredict}}, \code{\link{SignificantSubcluster}}
#' @export
#' @rdname GOpredict-methods
GOpredict<-function(Pfile, PPIdb, Gene2Annotation, p_value=0.001){
	#DataBase=read.delim(PPIdb, as.is=T)
	DataBase=PPIdb
	LenH=length(unique(rm.na(c(DataBase[,1],DataBase[,2]))))
	#goleafvectorlite=dget(Gene2Annotation)
	goleafvectorlite=Gene2Annotation
	orderup=Pfile
	proall=unique(c(orderup[,1], orderup[,2]))
	GOidall=vector(mode='list', length=length(proall))
	neibor =vector(mode='list', length=length(proall))
	aa=0
	for (i in proall){
		aa=aa+1
		t1=orderup[orderup[,1]==i,2]
		t2=orderup[orderup[,2]==i,1] 
		t=unique(rm.na(c(t1,t2)))
		GOid=character()
		for (k in t ){
			id= unique(goleafvectorlite[[k]])
			GOid= c(GOid,id)
		}
		GOidall[[aa]]=GOid
		neibor[[aa]]=length(t)
	}
	names(GOidall) <- proall
	names(neibor)  <- proall

	allgo=sort(table(unlist(goleafvectorlite)), decreasing=T)
	goname=names(allgo)
	gocount=unname(allgo)
	goLen=length(goleafvectorlite)
	allGanno=sum(gocount)
	GOannoAssign=character()
	ovlap=numeric()
	Pvalue=numeric()
	ind=0
	for (i in proall){
		ind=ind+1
		if (length(AA<-GOidall[[i]])>0){
			z=numeric()
			ov=numeric()
			for (k in unique(AA)){
				x1=c(gocount[match(k, goname)]-sum(AA==k), LenH-gocount[match(k, goname)]-(neibor[[i]]-sum(AA==k)))  # LenH is the number of human proteins in PPI data.
				x2=c(sum(AA==k), neibor[[i]]-sum(AA==k))
				xx=data.frame(x1=x1,x2=x2)
				z=c(z,fisher.test(xx)$p.value)
				ov=c(ov, sum(AA==k))
			}
			if (min(z)<p_value){
				GOannoAssign[ind]=unique(AA)[which.min(z)]
				Pvalue[ind]=min(z)
				ovlap[ind]=ov[which.min(z)]
			}
		}
	}

	TrueAnnoGO=rep(NA, length(proall))
	Ra=rep(NA, length(proall)) 
	Rao=rep(NA, length(proall)) 
	for (i in 1:length(GOannoAssign)){
		A<-GOannoAssign[i][!is.element(rm.na(GOannoAssign[i]),goleafvectorlite[[proall[i]]])]
		if (length(A)>0) {
			TrueAnnoGO[i]=A
			Ra[i]=paste(ovlap[i],'/', neibor[[i]], sep='')
			Rao[i]=ovlap[i]/neibor[[i]]
		}
	}

	X=which(!is.na(TrueAnnoGO))
	GOterm=character()
	for (i in TrueAnnoGO[X]){
		if (!is.null(GOTERM[[i]])){
			GOterm=c(GOterm, GOTERM[[i]]@Term)
		}
		else {GOterm=c(GOterm, GOSYNONYM[[i]]@Term)}
	}
	GOpredict=data.frame(Symbol=proall[X], GOID=TrueAnnoGO[X], GOterm=GOterm, Ratio=Ra[X], Pvalue=Pvalue[X])
	GOpredict
}
###GP=GOpredict(Pfile=OrderAll, PPIdb='SampleDatabase.txt', Gene2Annotation='GENE2GOtopLite', p_value=0.001)




### Section III: search for KEGG-term enriched subclusters (only for protein network) ###
#' @title Subclusters with KEGG annotations significantly enriched
#'
#' @description This function identifies subclusters whose members are significantly enrichment in certain KEGG pathways.
#'
#' @param Dendrogram An object in the class "dendrogram" returned from the function "ProteinCluster"
#' @param Gene2Annotation A list that maps KEGG pathway ID to genes. The names should be gene symbols and the elements should be KEGG pathways. i.e. 
#'
#'	$IGSF5
#'
#'	[1] "hsa04530" "hsa05120"
#'
#'	$OR2T8
#'
#'	[1] "hsa04740"
#'
#'	$hCG_1776980
#'
#'	[1] "hsa00020" "hsa00190" "hsa01100" "hsa05010" "hsa05012" "hsa05016"
#'
#'	......
#'
#' @param PPIdb A 2-column data frame consisting of binary interactions where each row i.e. c(A, B) represents an undirected edge (interaction) between gene A and gene B.  
#' @param KGremove If "TRUE", "hsa05200" and "hsa01100" will be excluded from KEGG-pathway based enrichment analysis as they are too broad.
#' @param SPoint The starting point for searching for the KEGG-pathway enriched subclusters
#' @param EPoint The endpoint for searching for the KEGG-pathway enriched subclusters
#' @param p_value Description for p_value
#' @param ini.p_value Description for ini.p_value
#' @usage SignificantSubcluster(Dendrogram, Gene2Annotation, PPIdb, KGremove=TRUE, 
#' SPoint=1, EPoint=9.7, p_value=0.001, ini.p_value=0.05)
#' @return This function identifies subclusters whose members are significantly enrichment in KEGG pathways and generate a dendrogram for these subclusters. 
#' @seealso \code{\link{SignificantPairs}}, \code{\link{ProteinCluster}}, \code{\link{KEGGpredict}}, \code{\link{GOpredict}}
#' @export
#' @rdname SignificantSubcluster-methods
SignificantSubcluster<-function(Dendrogram, Gene2Annotation, PPIdb, KGremove=TRUE, SPoint=1, EPoint=9.7, p_value=0.001, ini.p_value=0.05){

	#keggvectorR=dget(Gene2Annotation)
	keggvectorR=Gene2Annotation
	if (KGremove==TRUE){
		for (i in 1:length(keggvectorR)){       #first, remove 'hsa05200', 'hsa01100' from keggvector as they are too broad.
			a=match('hsa05200', keggvectorR[[i]])
			b=match('hsa01100', keggvectorR[[i]])
			c=rm.na(c(a,b))
			if (length(c)>0) {
				keggvectorR[[i]]=keggvectorR[[i]][-c]
			}
		}
	}
	allkegg=sort(table(unlist(keggvectorR)), decreasing=T)
	keggname=names(allkegg)
	keggcount=unname(allkegg)
	kgLen=length(keggvectorR)
	allKanno=sum(keggcount)
	
	#DataBase=read.delim(PPIdb, as.is=T)
	DataBase=PPIdb
	LenH=length(unique(rm.na(c(DataBase[,1],DataBase[,2]))))

	StartPoint=SPoint    # User should be able to pick from iniCut from 0.1-9.9 by themselves based on the cluster. 
	iniCut=StartPoint 
	cutden=cut(Dendrogram, iniCut)
	ini.P=list()
	ini.Kegg=list()
	ini.Pro=list()
	for (i in 1:length(cutden$lower)){
		if (length(order.dendrogram(cutden$lower[[i]]))>4) { 
			b=character()
			Pro=labels(cutden$lower[[i]])
			for (k in Pro){
				b=c(b, keggvectorR[[k]])
			}
			z=numeric()
			ov=numeric()
			if (length(b)>0) {
				for (kk in unique(b)){
					x1=c(keggcount[match(kk, keggname)]-sum(b==kk), LenH-keggcount[match(kk, keggname)]-(length(Pro)-sum(b==kk)))  #'keggcount' and 'keggname' are calculated in "I KEGG"; LenH is the number of human proteins in PPI data.
					x2=c(sum(b==kk), length(Pro)-sum(b==kk))
					xx=data.frame(x1=x1,x2=x2)
					z=c(z,fisher.test(xx)$p.value)
					ov=c(ov, sum(b==kk))
				}
				ini.P[[i]]=z
				ini.Kegg[[i]]=unique(b)
				ini.Pro[[i]]=Pro
			}
		}
	}

	ini.Thres=ini.p_value
	#ini.Thres=0.05/length(unlist(ini.P))
	pos.P=list()
	pos.Kegg=list()
	pos.Pro=list()
	s=0
	for (i in 1:length(ini.P)){
		if (length(ini.P[[i]])>0){
			if (any(ini.P[[i]]<ini.Thres)){
				s=s+1
				pos.P[[s]]=ini.P[[i]][which(ini.P[[i]]<ini.Thres)]
				pos.Kegg[[s]]=ini.Kegg[[i]][which(ini.P[[i]]<ini.Thres)]
				pos.Pro[[s]]=ini.Pro[[i]]
			}
		}
	}

	pri.P=pos.P
	pri.Kegg=pos.Kegg
	pri.Pro=pos.Pro
	ReP=pos.P
	EndPoint=EPoint
	CutOff=seq(iniCut+0.1, EndPoint, 0.1)
	oo=0

	for (u in CutOff){
		oo=oo+1
		cutden=cut(Dendrogram, u)
		Ind=1:length(cutden$lower)
		#p.value=numeric()
		#p.kegg=character()
		p.value=list()
		p.kegg=list()
		p.pro=list()
		for (i in 1:length(cutden$lower)){
			if (length(order.dendrogram(cutden$lower[[i]]))>4) {
				#print(i)
				b=character()
				Pro=labels(cutden$lower[[i]])
				for (k in Pro){
					b=c(b, keggvectorR[[k]])
				}
				z=numeric()
				ov=numeric()
				if (length(b)>0) {
					for (kk in unique(b)){
						x1=c(keggcount[match(kk, keggname)]-sum(b==kk), LenH-keggcount[match(kk, keggname)]-(length(Pro)-sum(b==kk)))  #'keggcount' and 'keggname' are calculated in "I KEGG"; LenH is the number of human proteins in PPI data.
						x2=c(sum(b==kk), length(Pro)-sum(b==kk))
						xx=data.frame(x1=x1,x2=x2)
						z=c(z,fisher.test(xx)$p.value)
						ov=c(ov, sum(b==kk))
					}
					#p.value[i]=z[which.min(z)]
					#p.kegg[i]=unique(b)[which.min(z)]
					p.value[[i]]=z
					p.kegg[[i]]=unique(b)
					p.pro[[i]]=Pro
				}
			}
		}
		names(p.value)<-1:length(p.value)
		names(p.kegg)<-1:length(p.value)
		names(p.pro)<-1:length(p.value)

		for (ii in 1:length(pos.Pro)){
			k1=unique(match.list(pos.Pro[[ii]], p.pro))
			k2=p.kegg[[k1]]
			k3=match(k2,pos.Kegg[[ii]])
			k4=pos.P[[ii]][rm.na(k3)]
			k5=p.value[[k1]][!is.na(k3)]
			if (any(k4>k5)){
				#print(ii)
				pos.P[[ii]][rm.na(k3)][which(k4>k5)]=k5[which(k4>k5)]
				ReP[[ii]][rm.na(k3)][which(k4>k5)]=u
			}
		}
		#print(oo)
	}

	for (i in 1:length(ReP)){
		ReP[[i]][which(ReP[[i]]<1)]=1
	}

	MaxID=integer()
	for (i in 1:length(pos.P)){
		MaxID[i]=which.min(pos.P[[i]])
	}
	Heit=integer()
	for (i in 1:length(ReP)){
		Heit[i]=ReP[[i]][MaxID[i]]
	}
	Kg=character()
	for (i in 1:length(pos.Kegg)){
		Kg[i]=pos.Kegg[[i]][MaxID[i]]
	}
	Pos=numeric()
	for (i in 1:length(pos.P)){
		Pos[i]=pos.P[[i]][MaxID[i]]
	}

	#make plots of subclusters.
	Ind=which(Pos< p_value)
	#pdf(SubGraphFile, height=GraphHit, width=GraphWid)
	SubLabel=list()
	a=0
	for (i in Ind){
		CD=cut(Dendrogram, Heit[i])
		PP=list()
		for (k in 1:length(CD$lower)){
			PP[[k]]=labels(CD$lower[[k]])
		}
		names(PP)=1:length(PP)
		PPP=CD$lower[[unique(match.list(pos.Pro[[i]], PP))]]
		a=a+1
		SubLabel[[a]]=labels(PPP)
		cex=min(30/length(labels(PPP)), 1.5)
		plot(PPP, nodePar=list(lab.cex=cex, cex=c(0,0), lab.col='tomato'), main=Kg[i], cex.main=3)
	}   # plot all subclusters that are significantly enriched in pathways.
	#dev.off()
}

#SignificantSubcluster(Dendrogram=dendMap, Gene2Annotation="GENE2KEGG", PPIdb=dfPPI, KGremove=TRUE, SPoint=1, EPoint=9.7)

