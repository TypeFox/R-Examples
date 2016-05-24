### R code from vignette source 'PANDA.Rnw'

###################################################
### code chunk number 1: PANDA.Rnw:34-36
###################################################
options(width=80)
options(continue=' ')


###################################################
### code chunk number 2: loadLibrary
###################################################
library(PANDA)


###################################################
### code chunk number 3: PANDA.Rnw:55-58
###################################################
data(GENE2GOtopLite)
data(GENE2KEGG)
data(KEGGID2NAME)
data(dfPPI)


###################################################
### code chunk number 4: PANDA.Rnw:63-65
###################################################
OrderAll=SignificantPairs(PPIdb=dfPPI)
head(OrderAll)


###################################################
### code chunk number 5: PANDA.Rnw:70-71
###################################################
dendMap=ProteinCluster(Pfile=OrderAll, Plot=TRUE, TextScaler=50)


###################################################
### code chunk number 6: PANDA.Rnw:76-78
###################################################
GP=GOpredict(Pfile=OrderAll, PPIdb=dfPPI, Gene2Annotation=GENE2GOtopLite, p_value=0.001)
head(GP)


###################################################
### code chunk number 7: PANDA.Rnw:81-83
###################################################
KP=KEGGpredict(Pfile=OrderAll, PPIdb=dfPPI, Gene2Annotation=GENE2KEGG, p_value=0.001, IDtoNAME=KEGGID2NAME)
KP


###################################################
### code chunk number 8: PANDA.Rnw:88-90
###################################################
SignificantSubcluster(Dendrogram=dendMap, Gene2Annotation=GENE2KEGG, 
	PPIdb=dfPPI, KGremove=TRUE, SPoint=1, EPoint=9.7)


###################################################
### code chunk number 9: sessionInfo
###################################################
sessionInfo()


