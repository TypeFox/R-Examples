### R code from vignette source 'ESEA.Rnw'

###################################################
### code chunk number 1: ESEA.Rnw:46-47
###################################################
library(ESEA)


###################################################
### code chunk number 2: ESEA.Rnw:58-64
###################################################
#obtain the data for background set of edges.
edgesbackgrand<-GetEdgesBackgrandData()
edgesbackgrand[1:10,]

#obtain the edge sets of pathways.
pathwayEdge.db<-GetPathwayEdgeData()


###################################################
### code chunk number 3: ESEA.Rnw:77-94
###################################################

#get example data
dataset<-GetExampleData("dataset")
class.labels<-GetExampleData("class.labels")
controlcharactor<-GetExampleData("controlcharactor")

#get the data for background set of edges
edgesbackgrand<-GetEdgesBackgrandData()

#Calculate the differential correlation score for edges
EdgeCorScore<-calEdgeCorScore(dataset, class.labels, controlcharactor, edgesbackgrand)

#print the top ten results to screen
EdgeCorScore[1:10]

#Each element is the differential correlation score of an edge and whose name
# correspond to the edge in the background set of edges.


###################################################
### code chunk number 4: ESEA.Rnw:104-151
###################################################
#get example data
dataset<-GetExampleData("dataset")
class.labels<-GetExampleData("class.labels")
controlcharactor<-GetExampleData("controlcharactor")

#get the data for background set of edges
edgesbackgrand<-GetEdgesBackgrandData()

#get the edge sets of pathways
pathwayEdge.db<-GetPathwayEdgeData()

#calculate the differential correlation score for edges
EdgeCorScore<-calEdgeCorScore(dataset, class.labels, controlcharactor, edgesbackgrand)

#identify dysregulated pathways by using the function ESEA.Main
Results<-ESEA.Main(
EdgeCorScore,
pathwayEdge.db,
weighted.score.type = 1, 
pathway = "kegg", 
gs.size.threshold.min = 15, 
gs.size.threshold.max = 1000,
reshuffling.type = "edge.labels",
nperm =10, 
p.val.threshold=-1,
FDR.threshold = 0.05, 
topgs =1
)

#print the summary results of pathways to screen
Results[[1]][[1]][1:5,]

#print the detail results of pathways to screen
Results[[2]][[1]][1:5,]

#write the summary results of pathways to tab delimited file.
write.table(Results[[1]][[1]], file = "kegg-SUMMARY RESULTS Gain-of-correlation.txt",
 quote=F, row.names=F, sep = "\t")
write.table(Results[[1]][[2]], file = "kegg-SUMMARY RESULTS Loss-of-correlation.txt",
 quote=F, row.names=F, sep = "\t")

#write the detail results of genes for each pathway with FDR.threshold< 0.05 to tab delimited file.
for(i in 1:length(Results[[2]])){
PathwayList<-Results[[2]][[i]]
filename <- paste(names(Results[[2]][i]),".txt", sep="", collapse="")
write.table(PathwayList, file = filename, quote=F, row.names=F, sep = "\t")
}


###################################################
### code chunk number 5: ESEA.Rnw:160-171
###################################################
#get example data
dataset<-GetExampleData("dataset")
class.labels<-GetExampleData("class.labels")
controlcharactor<-GetExampleData("controlcharactor")

#get the data for background set of edges
edgesbackgrand<-GetEdgesBackgrandData()

#calculate the differential correlation score for edges
EdgeCorScore<-calEdgeCorScore(dataset, class.labels, controlcharactor, edgesbackgrand)



###################################################
### code chunk number 6: GlobEdgeCorProfile
###################################################
#plot global edge correlation profile
PlotGlobEdgeCorProfile(EdgeCorScore)


###################################################
### code chunk number 7: ESEA.Rnw:191-212
###################################################
#get example data
dataset<-GetExampleData("dataset")
class.labels<-GetExampleData("class.labels")
controlcharactor<-GetExampleData("controlcharactor")

#get the data for background set of edges
edgesbackgrand<-GetEdgesBackgrandData()

#get the edge sets of pathways
pathwayEdge.db<-GetPathwayEdgeData()

#calculate the differential correlation score for edges
EdgeCorScore<-calEdgeCorScore(dataset, class.labels, controlcharactor,edgesbackgrand)

#identify dysregulated pathways by using the function ESEA.Main
#Results<-ESEA.Main(EdgeCorScore,pathwayEdge.db)
Results<-GetExampleData("PathwayResult")

#obtain the detail results of genes for a significant pathway
PathwayResult<-Results[[2]][1]



###################################################
### code chunk number 8: RunEdgeCorScore
###################################################
#Plot running edge enrichment score for the pathway result
PlotRunEnrichment(EdgeCorScore,PathwayResult,weighted.score.type = 1)


###################################################
### code chunk number 9: ESEA.Rnw:233-254
###################################################
#get example data
dataset<-GetExampleData("dataset")
class.labels<-GetExampleData("class.labels")
controlcharactor<-GetExampleData("controlcharactor")

#get the data for background set of edges
edgesbackgrand<-GetEdgesBackgrandData()

#get the edge sets of pathways
pathwayEdge.db<-GetPathwayEdgeData()

#calculate the differential correlation score for edges
EdgeCorScore<-calEdgeCorScore(dataset, class.labels, controlcharactor,edgesbackgrand)

#identify dysregulated pathways by using the function ESEA.Main
#Results<-ESEA.Main(EdgeCorScore,pathwayEdge.db)
Results<-GetExampleData("PathwayResult")

#obtain the detail results of genes for a significant pathway
PathwayNetwork<-Results[[2]][[1]]



###################################################
### code chunk number 10: PathwayNetwork
###################################################
#Plot the pathway-result network diagram, the edges which contribute to the ES are labeled with red.
PlotPathwayGraph(PathwayNetwork,layout=layout.random)


###################################################
### code chunk number 11: ESEA.Rnw:275-298
###################################################
#get example data
dataset<-GetExampleData("dataset")
class.labels<-GetExampleData("class.labels")
controlcharactor<-GetExampleData("controlcharactor")

#get the data for background set of edges
edgesbackgrand<-GetEdgesBackgrandData()

#get the edge sets of pathways
pathwayEdge.db<-GetPathwayEdgeData()

#calculate the differential correlation score for edges
EdgeCorScore<-calEdgeCorScore(dataset, class.labels, controlcharactor,edgesbackgrand)

#identify dysregulated pathways by using the function ESEA.Main
#Results<-ESEA.Main(EdgeCorScore,pathwayEdge.db)
Results<-GetExampleData("PathwayResult")

#obtain the detail results of genes for a significant pathway
PathwayNetwork<-Results[[2]][[1]]

#save the pathway-result network to a file which can be input to the Cytoscape software. 
SavePathway2File(PathwayNetwork,layout=layout.circle,file="Graph")


###################################################
### code chunk number 12: sessionInfo
###################################################
sessionInfo()


