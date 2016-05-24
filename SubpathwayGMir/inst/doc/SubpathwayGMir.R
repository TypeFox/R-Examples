### R code from vignette source 'SubpathwayGMir.Rnw'

###################################################
### code chunk number 1: SubpathwayGMir.Rnw:47-48
###################################################
library(SubpathwayGMir)


###################################################
### code chunk number 2: SubpathwayGMir.Rnw:55-60
###################################################
# get verified miRNA-target interactions
expMir2Tar <- GetK2riData("expMir2Tar")

# view first six rows of data
expMir2Tar[1:6,]


###################################################
### code chunk number 3: SubpathwayGMir.Rnw:71-82
###################################################
# get hsa-specificd miRNA-target interactions
expMir2Tar <- GetK2riData("expMir2Tar")
row1 <- which(expMir2Tar[["LowTHExps"]]=="YES")
row2 <- which(expMir2Tar[["Species"]]=="hsa")
relations <- unique(expMir2Tar[intersect(row1,row2),c(2:3)])

# get direct metabolic pathway graphs
 DirectGraphList <- GetK2riData("MetabolicGEGEEMGraph")

# get reconstructed direct pathway graph list
 DirectInteGraphList <- getInteGraphList(DirectGraphList, relations)  


###################################################
### code chunk number 4: DirectInteGraph
###################################################

# visualize the reconstructed direct pathway
plotGraph(DirectInteGraphList[[1]],layout=layout.random)


###################################################
### code chunk number 5: SubpathwayGMir.Rnw:106-110
###################################################
# get undirect metabolic pathway graphs
UndirectGraphList <- GetK2riData("MetabolicGEGEUEMGraph")
# get reconstructed undirect pathway graph list
UndirectInteGraphList <- getInteGraphList(UndirectGraphList, relations) 


###################################################
### code chunk number 6: UnDirectInteGraph
###################################################
# visualize the reconstructed undirect pathway
plotGraph(UndirectInteGraphList[[1]],layout=layout.random)


###################################################
### code chunk number 7: SubpathwayGMir.Rnw:139-146
###################################################
# get user-interested miRNAs and genes
moleculeList <- c(getBackground(type="gene")[1:1000],
               getBackground(type="miRNA")[1:2000])

# get located direct subpathways
DirectSubGraphList <- getLocSubGraph(moleculeList,DirectInteGraphList,
                      type="gene_miRNA",n=1,s=10)


###################################################
### code chunk number 8: DirectSubGraph
###################################################

# visualize the located direct pathway
plotGraph(DirectSubGraphList[[1]],layout=layout.random)


###################################################
### code chunk number 9: SubpathwayGMir.Rnw:172-175
###################################################
# get located undirect subpathways
UnDirectSubGraphList <- getLocSubGraph(moleculeList,UndirectInteGraphList,
                      type="gene_miRNA",n=1,s=10)


###################################################
### code chunk number 10: UnDirectSubGraph
###################################################

# visualize the located undirect pathway
plotGraph(UnDirectSubGraphList[[6]],layout=layout.random)


###################################################
### code chunk number 11: SubpathwayGMir.Rnw:205-211
###################################################
# identify significant direct subpathways
ann <- identifyGraph(moleculeList,DirectSubGraphList,type="gene_miRNA")
result <- printGraph(ann,detail=TRUE)

# view the result
head(result[,c(1:2,5:6)])


###################################################
### code chunk number 12: SubpathwayGMir.Rnw:217-226
###################################################
# identify significant undirect subpathways
ann <- identifyGraph(moleculeList,UnDirectSubGraphList,type="gene_miRNA")
result <- printGraph(ann,detail=TRUE)

# view the result
head(result[,c(1:2,5:6)])

# save the result
write.table(head(result),"result.txt",sep="\t",col.names=TRUE,row.names=FALSE)


###################################################
### code chunk number 13: SubpathwayGMir.Rnw:231-240
###################################################
# get verified miRNA-target interactions 
expMir2Tar <- GetK2riData(K2riData="expMir2Tar")

# get the background of miRNAs 
BGMiRNA <- GetK2riData(K2riData="BGMiRNA")

# get the background of genes 
BGGene <- GetK2riData(K2riData="BGGene")



###################################################
### code chunk number 14: SubpathwayGMir.Rnw:246-255
###################################################
# update the cel-specific environment variables 
 updateOrgEnvir("cel")

# show the current environment variables
 ls(k2ri)

# show the background of miRNAs
 k2ri$BGMiRNA[1:3]
 


###################################################
### code chunk number 15: sessionInfo
###################################################
sessionInfo()


