### R code from vignette source 'IntClustVignette.Rnw'

###################################################
### code chunk number 1: config
###################################################
options(width = 50)
options(continue=" ")
set.seed(1258794)


###################################################
### code chunk number 2: PackageAndData
###################################################
library(IntClust)
data(fingerprintMat)
data(targetMat)
data(geneMat)


###################################################
### code chunk number 3: SingleClust
###################################################
MCF7_F <- Cluster(Data=fingerprintMat,type="data",distmeasure="tanimoto", 
                  normalize=FALSE,method=NULL,clust="agnes",linkage="flexible",
                  gap=FALSE,maxK=55,StopRange=FALSE) 

MCF7_T <- Cluster(Data=targetMat,type="data",distmeasure="tanimoto",
                  normalize=FALSE,method=NULL,clust="agnes",linkage="flexible",
                  gap=FALSE,maxK=55,StopRange=FALSE)


###################################################
### code chunk number 4: SelectnrClusters
###################################################
List=list(fingerprintMat,targetMat)

NrClusters=SelectnrClusters(List=List,type="data",distmeasure=c("tanimoto",
                            "tanimoto"),nrclusters=seq(5,20),normalize=FALSE,
                             method=NULL,names=c("FP","TP"),StopRange=FALSE,
                             plottype="sweave",location=NULL)


###################################################
### code chunk number 5: Colors
###################################################
Colors <- ColorPalette(colors=c("chocolate","firebrick2", "darkgoldenrod2",
                       "darkgreen","blue2","darkorchid3","deeppink"),ncols=7) 


###################################################
### code chunk number 6: ClusterPlots
###################################################
ClusterPlot(Data1=MCF7_F,Data2=MCF7_F,nrcluster=7,cols=Colors,main="Clustering on 
            Fingerprints: Dendrogram",ylim=c(-0.1,1.8))

ClusterPlot(Data1=MCF7_T,Data2=MCF7_F,nrcluster=7,cols=Colors,main="Clustering on 
            Targets: Dendrogram",ylim=c(-0.1,2.5))


###################################################
### code chunk number 7: ClusterPlotF
###################################################
ClusterPlot(MCF7_F,MCF7_F,nrcluster=7,cols=Colors,main="Clustering on 
Fingerprints: Dendrogram",ylim=c(-0.1,1.8),plottype="sweave",location=NULL)


###################################################
### code chunk number 8: ClusterPlotT
###################################################
ClusterPlot(MCF7_T,MCF7_F,nrcluster=7,cols=Colors,main="Clustering on 
Targets: Dendrogram",ylim=c(-0.1,2.5),plottype="sweave",location=NULL) 


###################################################
### code chunk number 9: ADC
###################################################
L=list(fingerprintMat,targetMat)
MCF7_ADC=ADC(List=L,distmeasure="tanimoto",normalize=FALSE,method=NULL,
             clust="agnes", linkage="flexible")


###################################################
### code chunk number 10: ADEC
###################################################
L=list(fingerprintMat,targetMat)

MCF7_ADECa=ADECa(List=L,distmeasure="tanimoto",normalize=FALSE,method=NULL,
                 t=20,r=NULL,nrclusters=7,clust="agnes",linkage="flexible")
		 
MCF7_ADECb=ADECb(List=L,distmeasure="tanimoto",normalize=FALSE,method=NULL,
                 nrclusters=seq(5,25,1),clust="agnes",linkage="flexible")
		
MCF7_ADECc=ADECc(List=L,distmeasure="tanimoto",normalize=FALSE,method=NULL,
                 t=20,r=NULL, nrclusters=seq(5,25),clust="agnes",linkage="flexible")


###################################################
### code chunk number 11: CEC
###################################################
L=list(fingerprintMat,targetMat)

MCF7_CECa=CECa(List=L,distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,
               method=NULL,t=20,r=NULL,nrclusters=c(7,7),clust="agnes",
			   linkage=c("flexible","flexible"),StopRange=FALSE)
		 
MCF7_CECb=CECb(List=L,distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,
               method=NULL,nrclusters=seq(5,25,1),clust="agnes",
			   linkage=c("flexible","flexible"),StopRange=FALSE) 
		
MCF7_CECc=CECc(List=L,distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,
               method=NULL,t=20,r=NULL,nrclusters=seq(5,25),clust="agnes",
			   linkage=c("flexible","flexible"),StopRange=FALSE)		 


###################################################
### code chunk number 12: SNF
###################################################
L=list(fingerprintMat,targetMat)

MCF7_SNFa=SNFa(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,
               method=NULL,NN=5,mu=0.5,T=20,clust="agnes",linkage="ward",StopRange=FALSE)
		 
MCF7_SNFb=SNFb(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,
               method=NULL,NN=10,mu=0.5,T=20,clust="agnes",linkage="ward",StopRange=FALSE)
		
MCF7_SNFc=SNFc(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,
               method=NULL,NN=10,mu=0.5,T=20,clust="agnes",linkage="ward",StopRange=FALSE)	   


###################################################
### code chunk number 13: Weighted
###################################################
L=list(fingerprintMat,targetMat)

MCF7_Weighted=WeightedClust(L,type="data",distmeasure=c("tanimoto","tanimoto"),
                            normalize=FALSE,method=NULL,weight=seq(0,1,0.1),
                            WeightClust=0.5,clust="agnes",linkage="ward",
                            StopRange=FALSE)


###################################################
### code chunk number 14: WeightedSim
###################################################
L=list(fingerprintMat,targetMat)

MCF7_Weight=DetermineWeight_SimClust(List=L,type="data",weight=seq(1,0,by=-0.01),
                                     nrclusters=7,distmeasure=c("tanimoto","tanimoto"),
                                     normalize=FALSE,method=NULL,clust="agnes",
									 linkage=c("flexible","flexible"),gap=FALSE,maxK=50,
                                     names=c("FP","TP"),StopRange=FALSE,
                                     plottype="sweave",location=NULL)
					
L=list(MCF7_F,MCF7_T)	

MCF7_WeightedSim=WeightedSimClust(List=L,type="clusters",weight=MCF7_Weight$Weight,
                                  clust="agnes",linkage=c("flexible","flexible"),distmeasure=
                                  c("tanimoto","tanimoto"),normalize=FALSE,
                                  method=NULL,gap=FALSE,maxK=50,nrclusters=7,
                                  names=c("FP","TP"),AllClusters=FALSE,StopRange=FALSE)


###################################################
### code chunk number 15: WeightedSil (eval = FALSE)
###################################################
## L=list(fingerprintMat,targetMat)
## 
## MCF7_Weight=DetermineWeight_SilClust(List=L,type="data",weight=seq(0,1,by=0.01),
##                                      distmeasure=c("tanimoto","tanimoto"),
##                                      normalize=FALSE,nrclusters=7,names=c("FP","TP"),
##                                      nboot=1000,StopRange=FALSE,
##                                      plottype="sweave",location=NULL)


###################################################
### code chunk number 16: WonM
###################################################
L=list(fingerprintMat,targetMat)

MCF7_WonM=WonM(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),
               normalize=FALSE,method=NULL,nrclusters=seq(5,25),clust="agnes",
			   linkage=c("flexible","flexible"),StopRange=FALSE)


###################################################
### code chunk number 17: ComparePlotPrep1
###################################################
L=list(MCF7_F,MCF7_ADC,MCF7_ADECa,MCF7_ADECb,MCF7_ADECc,MCF7_CECa,MCF7_CECb,
       MCF7_CECc,MCF7_SNFa,MCF7_SNFb,MCF7_SNFc,MCF7_WonM,MCF7_Weighted,
       MCF7_WeightedSim,MCF7_T)

names=c("FP","ADC","ADECa","ADECb","ADECc","CECa","CECb","CECc","SNFa","SNFb",
        "SNFc","WonM","Weighted","WeightedSim","TP")


###################################################
### code chunk number 18: ComparePlot1
###################################################
ComparePlot(List=L,nrclusters=7,cols=Colors,fusionsLog=TRUE,WeightClust=TRUE,
            names=names,margins=c(9.1,4.1,4.1,4.1),plottype="sweave",location=NULL)


###################################################
### code chunk number 19: ComparePlotPrep2
###################################################
L=list(MCF7_F,MCF7_Weighted,MCF7_T)

names=c("FP",seq(1,0,-0.1),"TP")



###################################################
### code chunk number 20: ComparePlot2
###################################################
ComparePlot(List=L,nrclusters=7,cols=Colors,fusionsLog=TRUE,WeightClust=FALSE,
            names=names,margins=c(9.1,4.1,4.1,4.1),plottype="sweave",location=NULL)


###################################################
### code chunk number 21: CompareSvsMPrep (eval = FALSE)
###################################################
## LS=list(MCF7_F,MCF7_T)
## LM=list(MCF7_ADC,MCF7_ADECa,MCF7_ADECb,MCF7_ADECc,MCF7_CECa,MCF7_CECb,
## 		MCF7_CECc,MCF7_SNFa,MCF7_SNFb,MCF7_SNFc,MCF7_WonM,MCF7_Weighted,
## 		MCF7_WeightedSim)
## 
## nS=c("FP","TP")
## nM=c("ADC","ADECa","ADECb","ADECc","CECa","CECb","CECc","SNFa","SNFb",
##       "SNFc","WonM","Weighted","WeightedSim")
## CompareSvsM(ListS=LS,ListM=LM,nrclusters=7,cols=Colors,fusionsLogS=TRUE,
##             fusionsLogM=TRUE,WeightClustS=TRUE,WeightClustM=TRUE,
##             namesS=nS,namesM=nM,margins=c(8.1,3.1,3.1,4.1),plottype="sweave",
##             location=NULL)


###################################################
### code chunk number 22: CompareInter (eval = FALSE)
###################################################
## LS=list(MCF7_F,MCF7_T)
## LM=list(MCF7_ADC,MCF7_ADECa,MCF7_ADECb,MCF7_ADECc,MCF7_CECa,MCF7_CECb,
## 		MCF7_CECc,MCF7_SNFa,MCF7_SNFb,MCF7_SNFc,MCF7_WonM,MCF7_Weighted,
## 		MCF7_WeightedSim)
## 
## nS=c("FP","TP")
## nM=c("ADC","ADECa","ADECb","ADECc","CECa","CECb","CECc","SNFa","SNFb",
## 		"SNFc","WonM","Weighted","WeightedSim")
## Colors=c(Colors,"grey")
## CompareInteractive(ListM=LM,ListS=LS,nrclusters=7,cols=Colors,fusionsLogM=TRUE,
##                    fusionsLogS=TRUE,WeightClustM=TRUE,WeightClustS=TRUE,
##                    namesM=nM,namesS=nS,marginsM=c(8,4.5,2,2.5),marginsS
##                    =c(8,4.5,2,2.5),Interactive=TRUE,N=2)


###################################################
### code chunk number 23: CompareSilCluster (eval = FALSE)
###################################################
## List=list(fingerprintMat,targetMat)
## 
## Comparison=CompareSilCluster(List=List,type="data",distmeasure=c("tanimoto",
##                              "tanimoto"),normalize=FALSE,method=NULL,
## 		                      nrclusters=7,names=c("FP","TP"),nboot=1000,
##                               StopRange=FALSE,plottype="sweave",location=NULL)
## 					  
## Comparison


###################################################
### code chunk number 24: Comps
###################################################
L=list(MCF7_F,MCF7_ADC,MCF7_ADECa,MCF7_ADECb,MCF7_ADECc,MCF7_CECa,MCF7_CECb,
       MCF7_CECc,MCF7_SNFa,MCF7_SNFb,MCF7_SNFc,MCF7_WonM,MCF7_Weighted,
       MCF7_WeightedSim,MCF7_T)

Comps=FindCluster(List=L,nrclusters=7,select=c(15,4))
Comps


###################################################
### code chunk number 25: ClusterDistrPrep (eval = FALSE)
###################################################
## L=list(MCF7_F,MCF7_ADC,MCF7_ADECa,MCF7_ADECb,MCF7_ADECc,MCF7_CECa,MCF7_CECb,
##        MCF7_CECc,MCF7_SNFa,MCF7_SNFb,MCF7_SNFc,MCF7_WonM,MCF7_Weighted,
##        MCF7_WeightedSim,MCF7_T)
## 
## names=c("FP","ADC","ADECa","ADECb","ADECc","CECa","CECb","CECc","SNFa","SNFb",
## 		"SNFc","WonM","Weighted","WeightedSim","TP")
## 
## Tracking=TrackCluster(List=L,Selection=Comps,nrclusters=7,followMaxComps=FALSE,
##                       followClust=TRUE,fusionsLog=TRUE,WeightClust=TRUE,
##                       names=names,SelectionPlot=FALSE,Table=FALSE,
##                       CompleteSelectionPlot=TRUE,cols=Colors,plottype="sweave",
##                       location=NULL)


###################################################
### code chunk number 26: SimilarityHeatmap (eval = FALSE)
###################################################
## SimilarityHeatmap(Data=MCF7_Weighted$Clust,type="clust",distmeasure="tanimoto",
##                   normalize=FALSE,method="Q",cutoff=0.90,percentile=TRUE
##                   ,plottype="sweave",location=NULL)
## 	


###################################################
### code chunk number 27: HeatmapSelection (eval = FALSE)
###################################################
## CompsHeat=HeatmapSelection(Data=MCF7_Weighted$Clust$DistM,type="dist",
##                            cutoff=0.90,percentile=TRUE,dendrogram=
##                            MCF7_Weighted$Clust,width=7,height=7)
## 
## CompsHeat


###################################################
### code chunk number 28: DiffGenes
###################################################
MCF7_Genes=DiffGenes(List=NULL,Selection=Comps,GeneExpr=geneMat,
		             nrclusters=7,method="limma",sign=0.05,
                     topG=10,fusionsLog=TRUE,WeightClust=TRUE,names=NULL)
Genes=MCF7_Genes$Selection$Genes$TopDE$ID				  


###################################################
### code chunk number 29: ProfilePlot
###################################################
ProfilePlot(Genes=Genes[1:5],Comps=Comps,GeneExpr=geneMat,
            Raw=FALSE,Order=MCF7_F,Color=MCF7_F,nrclusters=7,
            Clusters=NULL,cols=Colors,AddLegend=TRUE,
            margins=c(8.1,4.1,1.1,6.5),cex=0.75,plottype="sweave",
            location=NULL)



###################################################
### code chunk number 30: Pathways (eval = FALSE)
###################################################
## library(MLP)
## data(GeneInfo)
## data(GS)
## L=list(MCF7_Genes)
## MCF7_Paths=PathwaysIter(List=L,Selection=Comps,names=NULL,
##                         GeneExpr=geneMat,nrclusters=7,method=
##                         c("limma", "MLP"),GeneInfo=GeneInfo,
##                         geneSetSource = "GOBP",topP=NULL,topG=NULL,
##                         GENESET=GS,sign=0.05,niter=2,
##                         fusionsLog=TRUE,WeightClust=TRUE)
## 
## MCF7_Paths_Inter=Geneset.intersect(MCF7_Paths,Selection=TRUE,sign=0.05,
##                                    names=NULL,seperatetables=FALSE,
##                                    separatepvals=FALSE)
## 


###################################################
### code chunk number 31: ChooseFeatures
###################################################
MCF7_Feat=ChooseCluster(Interactive=FALSE,LeadCpds=Comps,ClusterResult=MCF7_F,
                        ColorLab=MCF7_F,BinData=list(fingerprintMat,
                        targetMat),ContData=NULL,Datanames=c("FP","TP"),geneMat,topChar = 20,
                        topG = 20,sign=0.05,nrclusters=7,N=1)
				

					


###################################################
### code chunk number 32: FeaturesPlotFP (eval = FALSE)
###################################################
## BinFeaturesPlot(LeadCpds=Comps,OrderLab=MCF7_F,Features=
##                 MCF7_Feat$Characteristics$FP$TopFeat$Names,Data=fingerprintMat,
##                 ColorLab=MCF7_F,nrclusters=7,cols=Colors,name=c("FP"),plottype="sweave",
##                 location=NULL)		


###################################################
### code chunk number 33: FeaturesPlotT (eval = FALSE)
###################################################
## BinFeaturesPlot(LeadCpds=Comps,OrderLab=MCF7_F,Features=
##                 MCF7_Feat$Characteristics$TP$TopFeat$Names,Data=targetMat,ColorLab=MCF7_F,
##                 nrclusters=7,cols=Colors,name=c("TP"),plottype="sweave",
##                 location=NULL)


###################################################
### code chunk number 34: FeaturesPlotFP1
###################################################
BinFeaturesPlot(LeadCpds=Comps,OrderLab=MCF7_F,Features=
                MCF7_Feat$Characteristics$FP$TopFeat$Names,Data=fingerprintMat,
                ColorLab=MCF7_F,nrclusters=7,cols=Colors,name=c("FP"),
                plottype="sweave",location=NULL)


###################################################
### code chunk number 35: FeaturesPlotT1
###################################################
BinFeaturesPlot(LeadCpds=Comps,OrderLab=MCF7_F,Features=
                MCF7_Feat$Characteristics$TP$TopFeat$Names,Data=targetMat,
                ColorLab=MCF7_F,nrclusters=7,cols=Colors,name=c("TP"),
                plottype="sweave",location=NULL)


###################################################
### code chunk number 36: sessionInfo
###################################################
toLatex(sessionInfo())


