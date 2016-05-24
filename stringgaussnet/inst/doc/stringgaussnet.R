## ----eval=FALSE----------------------------------------------------------
#  library(stringgaussnet)
#  
#  data(SpADataExpression) # Import example expression data
#  data(SpADEGenes) # Import example DE genes analysis results
#  data(SpASamples) # Import example sample description
#  # We firstly import all example data from stringgaussnet package.
#  # SpASamples is not compulsory for using stringgaussnet, but is useful for creating a
#  # factor for subsetting gaussian networks generation.
#  
#  head(SpASamples,5)

## ----eval=FALSE----------------------------------------------------------
#  # We can see here what sample descriptions look like. LPStime is the LPS stimulation
#  # duration for MD-DCs.
#  SpAData<-DEGeneExpr(t(SpADataExpression),SpADEGenes)
#  
#  print(SpAData,5)

## ----eval=FALSE----------------------------------------------------------
#  # Here we see that we have 57 samples and 75 DE genes.
#  
#  # There is a preview of the expression data. Row names correspond to the column chipnum
#  # in SpAsamples.
#  
#  # For each gene in row, we have attributes that will be added further as node attributes
#  # in our networks. Notably, we have here the both gene identifiers HGNC symbols and
#  # Ensembl IDs, and p-values and fold changes computed by LIMMA.

## ----eval=FALSE----------------------------------------------------------
#  library(stringgaussnet)
#  data(SpADataExpression)
#  data(SpADEGenes)
#  SpAData<-DEGeneExpr(t(SpADataExpression),SpADEGenes)
#  SpASTRINGNet<-getSTRINGNet(SpAData)
#  # Here we get the STRING network with default parameters.
#  # You can type help(“getSTRINGNet”) for more details.
#  
#  print(SpASTRINGNet,5)

## ----eval=FALSE----------------------------------------------------------
#  # We can see that STRING gave a network with 53 from the 75 initial genes
#  # entered in the API. 183 additional nodes were used to construct the network.
#  
#  # We can see also that 17 interactions between initial nodes were found, on contrary to
#  # 5168 including added nodes. Multiple edges are not taken into account, which means
#  # that the print function displays unique pairs of genes as interactions.
#  
#  # We have here, for each edge, an interaction source attribute and the corresponding
#  # score given by STRING. Combined scores are also entered, with the label “combined”.
#  
#  # As we can see, differential analysis results are used as node attributes.
#  
#  summary(SpASTRINGNet)

## ----eval=FALSE----------------------------------------------------------
#  # Here we have score summaries for each interaction source and by making a distinction
#  # between initial and added nodes.
#  
#  PPISpASTRINGNet <- selectInteractionTypes(SpASTRINGNet,
#    c("coexpression","experimental","knowledge"), 0.9)
#  # Here we select only interactions of kind “coexpression”, “experimental” and
#  # “knowledge”, with a score filtering threshold of 0.9.
#  
#  print(PPISpASTRINGNet,5)

## ----eval=FALSE----------------------------------------------------------
#  summary(PPISpASTRINGNet)

## ----eval=FALSE----------------------------------------------------------
#  library(stringgaussnet)
#  data(SpADataExpression)
#  data(SpADEGenes)
#  SpAData<-DEGeneExpr(t(SpADataExpression),SpADEGenes)
#  SpASTRINGNet<-getSTRINGNet(SpAData)
#  PPISpASTRINGNet <- selectInteractionTypes(SpASTRINGNet,
#    c("coexpression","experimental","knowledge"), 0.9)
#  shortPathSpANet<-getShortestPaths(PPISpASTRINGNet)
#  # Here we get the short paths STRING network with default parameters.
#  # You can type help(“getShortestPaths”) for more details.
#  shortPathSpANet<-FilterEdges(shortPathSpANet,5)
#  # Here we can filter on the distance between two nodes.
#  print(shortPathSpANet,5)

## ----eval=FALSE----------------------------------------------------------
#  # Here we don't have any added node, because we summarized the network only between
#  # initial nodes. We can see that 20 genes were used as intermediates.
#  
#  # We have unique edges with distance, number of intermediates and intermediate names
#  # as attributes.

## ------------------------------------------------------------------------
library(stringgaussnet)
data(SpADataExpression)
data(SpADEGenes)
SpAData<-DEGeneExpr(t(SpADataExpression),SpADEGenes)
NodesForSIMoNe<-rownames(SpADEGenes)[1:17]
GaussianSpAData<-DEGeneExpr(
  t(SpADataExpression[NodesForSIMoNe,]),
  SpADEGenes[NodesForSIMoNe,])
# We select a reasonable number of genes for SIMoNe network inference.
# We advice to take a number of genes being inferior to the sample size.
#pickSIMoNeParam(GaussianSpAData)
# We use a series of plot provided with the simone package to see which penalty level
# we can use for the graphical LASSO regression.
GlobalSIMoNeNet<-getSIMoNeNet(GaussianSpAData)
# Here we get the SIMoNe network with default parameters.
# You can type help(“getSIMoNeNet”) for more details.
GlobalSIMoNeNet<-FilterEdges(GlobalSIMoNeNet,0.4)
# Here we can filter on the absolute values of rho being superior to 0.4.
print(GlobalSIMoNeNet,5)
# We have unique edges with Theta scores, and spearman's rhos and p-values.
plot(GlobalSIMoNeNet)
# Here we have the network displayed with the plot function provided with the
# simone package.

## ------------------------------------------------------------------------
library(stringgaussnet)
data(SpADataExpression)
data(SpADEGenes)
SpAData<-DEGeneExpr(t(SpADataExpression),SpADEGenes)
NodesForSIMoNe<-rownames(SpADEGenes)[1:17]
GaussianSpAData<-DEGeneExpr(
  t(SpADataExpression[NodesForSIMoNe,]),
  SpADEGenes[NodesForSIMoNe,])
#pickWGCNAParam(GaussianSpAData)
# Here we use a list of plots to help in choosing the right parameter for
# WGCNA computing. You can see help(“pickWGCNAParam”) for more details.
GlobalWGCNANet<-getWGCNANet(GaussianSpAData)
# Here we get the WGCNA network with default parameters.
# You can type help(“getWGCNANet”) for more details.
print(GlobalWGCNANet,5)
# We have adjacency scores, and spearman's rhos and p-values for each edge.
plot(GlobalWGCNANet)
# Here we have the network displayed with the plot function provided in the
# simone package.

## ----eval=FALSE----------------------------------------------------------
#  library(stringgaussnet)
#  data(SpADataExpression)
#  data(SpADEGenes)
#  SpAData<-DEGeneExpr(t(SpADataExpression),SpADEGenes)
#  NodesForSIMoNe<-rownames(SpADEGenes)[1:17]
#  GaussianSpAData<-DEGeneExpr(
#    t(SpADataExpression[NodesForSIMoNe,]),
#    SpADEGenes[NodesForSIMoNe,])
#  GlobalSIMoNeNet<-getSIMoNeNet(GaussianSpAData,
#    AddAnnotations=TRUE)
#  # Here we use the parameter “AddAnnotations=TRUE” to add annotations to genes from
#  # the network.
#  print(GlobalSIMoNeNet,5)

## ----eval=FALSE----------------------------------------------------------
#  # We can see that we have gene annotations added by biomaRt and gene product
#  # localizations provided by Gene Ontology.

## ------------------------------------------------------------------------
library(stringgaussnet)
data(SpADataExpression)
data(SpADEGenes)
data(SpASamples)
SpAData<-DEGeneExpr(t(SpADataExpression),SpADEGenes)
StatusFactor<-SpASamples$status
names(StatusFactor)<-SpASamples$chipnum
# We create a factor vector based on the status.
NodesForSIMoNe<-rownames(SpADEGenes)[1:17]
GaussianSpAData<-DEGeneExpr(
  t(SpADataExpression[NodesForSIMoNe,]),
  SpADEGenes[NodesForSIMoNe,])
StatusFactorSIMoNeNet<-FactorNetworks(GaussianSpAData,
  StatusFactor,"SIMoNe")
# We infer different SIMoNe networks on different groups of samples
# (patients and controls).
StatusFactorSIMoNeNet<-FilterEdges(StatusFactorSIMoNeNet,0.4)
# We can filter on edges, like for SIMoNeNet.
print(StatusFactorSIMoNeNet)
# We can see that we have a preview of inferred networks for each level.
par(mfrow=c(2,1))
plot(StatusFactorSIMoNeNet,interactiveMode=F)
# Here we can display networks inferred in patients and controls specifically, like in
# the Figure 5, with the provided function in the simone package.

#compareFactorNetworks(StatusFactorSIMoNeNet)
# Here we can have a series of plots to compare results of inferred networks for
# each level. You can see help (“compareFactorNetworks”) for more details.

#StatusFactorSIMoNeNet<-FactorNetworks(GaussianSpAData,
#  StatusFactor,"SIMoNe",list(AddAnnotations=TRUE))
# This is how you should do if you wanted to add gene annotations.

## ----eval=FALSE----------------------------------------------------------
#  library(stringgaussnet)
#  data(SpADataExpression)
#  data(SpADEGenes)
#  data(SpASamples)
#  SpAData<-DEGeneExpr(t(SpADataExpression),SpADEGenes)
#  StatusFactor<-SpASamples$status
#  names(StatusFactor)<-SpASamples$chipnum
#  NodesForSIMoNe<-rownames(SpADEGenes)[1:17]
#  GaussianSpAData<-DEGeneExpr(
#    t(SpADataExpression[NodesForSIMoNe,]),
#    SpADEGenes[NodesForSIMoNe,])
#  MultiSpAData<-MultiDEGeneExpr(GaussianSpAData,
#    DEGeneExpr(t(SpADataExpression[18:34,]),
#    SpADEGenes[18:34,]),
#    DEGeneExpr(t(SpADataExpression[35:51,]),
#    SpADEGenes[35:51,]))
#  # We create multiple lists of DE genes results, and then a list of DEGeneExpr objects,
#  # by subsetting the original data.
#  print(MultiSpAData)

## ----eval=FALSE----------------------------------------------------------
#  # We have, by the specific print function, a preview of each DEGeneExpr object
#  # in the list.
#  MultiSpANetworks<-MultiNetworks(MultiSpAData,
#    SelectInteractionsSTRING=c("coexpression",
#    "experimental","knowledge"),STRINGThreshold=0.9,
#    FilterSIMoNeOptions=list(Threshold=0.4),
#    Factors=StatusFactor)
#  # We create an object of class MultiNetworks, which allows to generate all kinds of
#  # networks in multiple lists of DEGeneExpr objects in one line.
#  print(MultiSpANetworks)

## ----eval=FALSE----------------------------------------------------------
#  # We simply have a summary of the used method to generate the MultiNetworks object.
#  
#  #MultiSpANetworks<-MultiNetworks(MultiSpAData,
#  #  SelectInteractionsSTRING=c("coexpression",
#  #  "experimental","knowledge"),STRINGThreshold=0.9,
#  #  FilterSIMoNeOptions=list(Threshold=0.4),
#  #  Factors=StatusFactor,
#  #  STRINGOptions=list(AddAnnotations=TRUE),
#  #  SIMoNeOptions=list(AddAnnotations=TRUE),
#  #  WGCNAOptions=list(AddAnnotations=TRUE))
#  # This is how you should do if you wanted to add gene annotations for all network generations.

## ----eval=FALSE----------------------------------------------------------
#  library(stringgaussnet)
#  data(SpADataExpression)
#  data(SpADEGenes)
#  SpAData<-DEGeneExpr(t(SpADataExpression),SpADEGenes)
#  NodesForSIMoNe<-rownames(SpADEGenes)[1:17]
#  GaussianSpAData<-DEGeneExpr(
#    t(SpADataExpression[NodesForSIMoNe,]),
#    SpADEGenes[NodesForSIMoNe,])
#  GlobalSIMoNeNet<-getSIMoNeNet(GaussianSpAData,
#    AddAnnotations=TRUE)
#  GlobalSIMoNeNet<-FilterEdges(GlobalSIMoNeNet,0.4)
#  export(GlobalSIMoNeNet,YourDirPath)
#  # Replace YourDirPath by the directory where will be saved edge and node attributes.
#  # If the directory exists, this will not be overwritten,
#  # excepted if you use the parameter “overwrite=TRUE”.

## ----eval=FALSE----------------------------------------------------------
#  library(stringgaussnet)
#  data(SpADataExpression)
#  data(SpADEGenes)
#  data(SpASamples)
#  SpAData<-DEGeneExpr(t(SpADataExpression),SpADEGenes)
#  StatusFactor<-SpASamples$status
#  names(StatusFactor)<-SpASamples$chipnum
#  NodesForSIMoNe<-rownames(SpADEGenes)[1:17]
#  GaussianSpAData<-DEGeneExpr(
#    t(SpADataExpression[NodesForSIMoNe,]),
#    SpADEGenes[NodesForSIMoNe,])
#  MultiSpAData<-MultiDEGeneExpr(GaussianSpAData,
#    DEGeneExpr(t(SpADataExpression[18:34,]),
#    SpADEGenes[18:34,]),
#    DEGeneExpr(t(SpADataExpression[35:51,]),
#    SpADEGenes[35:51,]))
#  MultiSpANetworks<-MultiNetworks(MultiSpAData,
#    SelectInteractionsSTRING=c("coexpression",
#    "experimental","knowledge"),STRINGThreshold=0.9,
#    FilterSIMoNeOptions=list(Threshold=0.4),
#    Factors=StatusFactor,
#    STRINGOptions=list(AddAnnotations=TRUE),
#    SIMoNeOptions=list(AddAnnotations=TRUE),
#    WGCNAOptions=list(AddAnnotations=TRUE))
#  resetCytoscapeSession()
#  # We reset and create an empty session in Cytoscape.
#  # Please be sure that Cytoscape is running with the cyREST plugin installed.
#  addMultiGraphToCytoscape(MultiSpANetworks,
#    points.size.map="P.Value",points.fill.map="logFC")
#  # We add automatically all generated networks in the Cytoscape session,
#  # with predefined styles.
#  saveCytoscapeSession(YourFilePath)
#  # We can save the current Cytoscape session in a .cys file.
#  # Replace YourFilePath with the path where you would like to save.
#  # The .cys extension is automatically added if necessary.

