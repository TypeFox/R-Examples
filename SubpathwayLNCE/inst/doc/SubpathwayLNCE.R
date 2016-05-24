### R code from vignette source 'SubpathwayLNCE.Rnw'

###################################################
### code chunk number 1: SubpathwayLNCE.Rnw:47-48
###################################################
library(SubpathwayLNCE)


###################################################
### code chunk number 2: SubpathwayLNCE.Rnw:55-64
###################################################
#obtain the data for candidate lncRNA-mRNA interaction.
interaction<-GetExampleData(exampleData="pp")
# view first six rows of data
interaction[1:6,]
#obtain the data for undirect KEGG metabolic pathway graphs with genes as nodes
 g2<-GetExampleData(exampleData="g2")
#obtain example data of mathed mRNA-lncRNA expression profiles
 #GeneExp<-GetExampleData(exampleData="GeneExp")
 #LncExp<-GetExampleData(exampleData="LncExp")


###################################################
### code chunk number 3: SubpathwayLNCE.Rnw:75-85
###################################################
#obtain example data of mathed mRNA-lncRNA expression profiles
 GeneExp<-GetExampleData(exampleData="GeneExp")
 LncExp<-GetExampleData(exampleData="LncExp")
#calculated co-expression coefficient,the significant positive threshold is 0.025
 LncGenePairs<-getLncGenePairs(GeneExp,LncExp,a=0.025)
#obtain the data for undirect KEGG metabolic pathway graphs with genes as nodes
 g2<-GetExampleData(exampleData="g2")
# get reconstructed undirect pathway graph list
 #interUMGraph<-getInteGraphList(g2,LncGenePairs)
 


###################################################
### code chunk number 4: SubpathwayLNCE.Rnw:92-106
###################################################
#obtain the data for undirect KEGG metabolic pathway graphs with genes as nodes
  g2<-GetExampleData(exampleData="g2")
 #obtain example data of mathed mRNA-lncRNA expression profiles
 #GeneExp<-GetExampleData(exampleData="GeneExp")
 #LncExp<-GetExampleData(exampleData="LncExp")
#calculated co-expression coefficient,the significant positive threshold is 0.025
 #LncGenePairs<-getLncGenePairs(GeneExp,LncExp,a=0.025)
# get reconstructed undirect pathway graph list
# To improve efficiency, a fraction of signal pathway as case
  LncGenePairs<-GetExampleData(exampleData="LncGenePairs")
  interUMGraph<-getInteGraphList(g2[42:45],LncGenePairs)
### Integrate lncRNAs of competitive regulation into KEGG pathway graphs ###
  ##LncGenePairs<-GetExampleData(exampleData="LncGenePairs")
  ##inteUMGraph<-getInteUMGraph(LncGenePairs)


###################################################
### code chunk number 5: UnDirectInteGraph
###################################################
# visualize the reconstructed undirect pathway
  #LncGenePairs<-GetExampleData(exampleData="LncGenePairs")
  #inteUMGraph<-getInteUMGraph(LncGenePairs)
  plotGraphL(interUMGraph[[1]],vertex.label=getNodeLabel)


###################################################
### code chunk number 6: SubpathwayLNCE.Rnw:132-144
###################################################
### Integrate lncRNAs of competitive regulation into KEGG pathway graphs ###
  LncGenePairs<-GetExampleData(exampleData="LncGenePairs")
  #inteUMGraph<-getInteUMGraph(LncGenePairs)
  # To improve efficiency, a fraction of signal pathway as case
  LncGenePairs<-GetExampleData(exampleData="LncGenePairs")
  interUMGraph<-getInteGraphList(g2[42:45],LncGenePairs)
### get user-interested lncRNAs and genes sets.
##geneLnc<-c(getBackground(type="gene")[1:3000],unique(LncGenePairs[1,]))
  geneLnc<-GetExampleData(exampleData="geneLnc")
# get locate subpathways.
  sub<-getLocSubGraphLnc(geneLnc,interUMGraph,type="gene_lncRNA",n=1,s=8)
  


###################################################
### code chunk number 7: SubpathwayLNCE.Rnw:152-166
###################################################
### Integrate lncRNAs of competitive regulation into KEGG pathway graphs ###
  #LncGenePairs<-GetExampleData(exampleData="LncGenePairs")
  #inteUMGraph<-getInteUMGraph(LncGenePairs)
### get user-interested lncRNAs and genes sets.
##geneLnc<-c(getBackground(type="gene")[1:3000],unique(LncGenePairs[1,]))
  geneLnc<-GetExampleData(exampleData="geneLnc")
# get locate subpathways.
  #sub<-getLocSubGraphLnc(geneLnc,interUMGraph,type="gene_lncRNA",n=1,s=8)
  sub<-GetExampleData(exampleData="sub")
  # To improve efficiency, a fraction of signal subpathway as case
  SubcodeLncResult<-identifyLncGraphW(geneLnc,sub[50:55],type="gene_lncRNA",bet=1)
  #SubcodeLncResult<-identifyLncGraphW(geneLnc,sub,type="gene_lncRNA",bet=1)
  #resultT<-printGraphW(SubcodeLncResult,detail=TRUE)
  #write.table(resultT,file="result.txt",sep="\t",row.names=F,quote=F)


###################################################
### code chunk number 8: PlotAnnGraph
###################################################
  plotAnnGraph("path:04916_1",sub,SubcodeLncResult,gotoKEGG=FALSE,vertex.label=getNodeLabel)


###################################################
### code chunk number 9: sessionInfo
###################################################
sessionInfo()


