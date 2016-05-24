########################################################################
# Categorical Network Class

setClass("catNetwork",  
         representation(
                        objectName = "character",
                        numnodes   = "integer", 
                        nodes      = "vector",
                        meta       = "character",
                        maxpars    = "integer", pars="list",  
                        maxcats    = "integer", cats="list",
                        probs      = "list",
                        complx     = "integer",
                        nodecomplx = "ANY",
                        loglik     = "numeric",
                        nodelik    = "ANY",
                        nodeSampleSizes = "vector"), 
         validity = function(object) validCatNetwork(object),
         package = "sdnet"
	   )

setGeneric("cnAddNode", 
          function(object, node, cats, pars, prob)
           standardGeneric("cnAddNode")
           )

setGeneric("cnNumNodes", function(object)
           standardGeneric("cnNumNodes")
           )

# returns list of nodes
setGeneric("cnNodes", function(object, which)
           standardGeneric("cnNodes")
           )

## returns list of edges
setGeneric("cnEdges", function(object, which)
           standardGeneric("cnEdges")
           )

setGeneric("cnMatEdges", function(object)
           standardGeneric("cnMatEdges")
           )

## returns pars sets as a list
setGeneric("cnParents", function(object, which)
           standardGeneric("cnParents")
           )

## returns pars sets as a matrix
setGeneric("cnMatParents", 
          function(object, nodeorder)
           standardGeneric("cnMatParents")
           )

setGeneric("cnProb", function(object, which=NULL)
           standardGeneric("cnProb")
           )

setGeneric("cnPlotProb", function(object, which=NULL)
           standardGeneric("cnPlotProb")
           )

setGeneric("cnDot", function(object, file=NULL, format="ps", nodestyle=NULL, edgestyle=NULL)
           standardGeneric("cnDot")
           )

setGeneric("cnSamples", function(object, numsamples=1, pert=NULL, output="frame", as.index=FALSE, naRate=0)
           standardGeneric("cnSamples")
           )

setGeneric("cnSetProb", function(object, data, pert=NULL, nodeCats=NULL, softmode=FALSE)
           standardGeneric("cnSetProb")
           )

setGeneric("cnLoglik", function(object, data, pert=NULL, bysample=FALSE, softmode=FALSE)
           standardGeneric("cnLoglik")
           )

setGeneric("cnNodeLoglik", function(object, node, data, pert=NULL, softmode=FALSE, klmode=FALSE)
           standardGeneric("cnNodeLoglik")
           )

setGeneric("cnNodeLoglikError", function(object, node, data, pert=NULL)
           standardGeneric("cnNodeLoglikError")
           )

setGeneric("cnPredict", function(object, data)
           standardGeneric("cnPredict")
           )
           
setGeneric("cnEvaluate", function(object, data, pert=NULL, maxParentSet = 0, maxcomplx=0, echo=FALSE)
           standardGeneric("cnEvaluate")
           )

setGeneric("cnComplexity", function(object, node=NULL, include.unif=TRUE)
           standardGeneric("cnComplexity")
           )

setGeneric("cnKLComplexity", function(object, node=NULL)
           standardGeneric("cnKLComplexity")
           )

## returns a graph object
setGeneric("as.graph", function(object)
           standardGeneric("as.graph")
           )

setGeneric("cnJointProb", 
          function(object,nodes)
           standardGeneric("cnJointProb")
           )

setGeneric("cnCondProb", 
          function(object,x,y)
           standardGeneric("cnCondProb")
           )

setGeneric("cnJointKLdist", 
          function(object1, object2, ...)
           standardGeneric("cnJointKLdist")
           )

setGeneric("cnNodeMarginalProb", 
          function(object, node)
           standardGeneric("cnNodeMarginalProb")
           )

setGeneric("cnMarginalKLdist", 
          function(object1, object2,...)
           standardGeneric("cnMarginalKLdist")
           )

setGeneric("cnMarginalKLdistList", 
          function(object1, object2list,...)
           standardGeneric("cnMarginalKLdistList")
           )

setGeneric("cnCondKLdist", 
          function(object1, object2, ...)
           standardGeneric("cnCondKLdist")
           )

setGeneric("cnCompare", 
          function(object1, object2, extended = FALSE)
           standardGeneric("cnCompare")
           )

setGeneric("cnSubNetwork", 
          function(object, nodeIndices, indirectEdges = FALSE)
           standardGeneric("cnSubNetwork")
           )

setGeneric("cnReorderNodes", 
          function(object, nodeIndices)
           standardGeneric("cnReorderNodes")
           )

setGeneric("cnOrder", 
          function(object)
           standardGeneric("cnOrder")
           )

setGeneric("cnCluster", 
          function(object)
           standardGeneric("cnCluster")
           )

setGeneric("cnClusterSep", 
          function(object, data, pert=NULL)
           standardGeneric("cnClusterSep")
           )

# draw a graph
setGeneric("cnPlot", 
          function(object, file=NULL)
           standardGeneric("cnPlot")
           )

setGeneric("cnFindAIC", function(object, numsamples)
           standardGeneric("cnFindAIC")
           )

setGeneric("cnFindBIC", function(object, numsamples)
           standardGeneric("cnFindBIC")
           )

setGeneric("cnFindKL", function(object, numsamples)
           standardGeneric("cnFindKL")
           )

setGeneric("addNetworkNode", 
          function(object, newnode, newparents, loglik, problist, nodecats, ...)
           standardGeneric("addNetworkNode")
           )

setGeneric("replaceNetworkNode", 
          function(object, newnode, newparents, loglik, problist, ...)
           standardGeneric("replaceNetworkNode")
           )

setGeneric("updateNetworkNode", 
          function(newnet, listnet, ...)
           standardGeneric("updateNetworkNode")
           )

setGeneric("cnPearsonTest", 
          function(object, data)
           standardGeneric("cnPearsonTest")
           )

setGeneric("cnMarParents", 
          function(object, flags = NULL)
           standardGeneric("cnMarParents")
           )

#########################################################################
## Categorical Network Diagnostic Class

setClass("catNetworkDistance",  
         representation(
                        hamm    = "numeric",
                        hammexp = "ANY",
                        tp = "numeric",
                        fp = "numeric",
                        fn = "numeric",
                        pr = "numeric",
                        sp = "numeric",
                        sn = "numeric",
                        fscore    = "numeric",
                        skel.tp   = "numeric",
                        skel.fp   = "numeric",
                        skel.fn   = "numeric",
                        order.fp  = "ANY",
                        order.fn  = "ANY", 
                        markov.fp = "ANY",
                        markov.fn = "ANY",
                        KLdist    = "ANY"),
         package = "sdnet"
	   )

setClass("catNetworkEvaluate",  
         representation(
                        numnodes   = "integer",
                        numsamples = "integer", 
                        nets    = "list",
                        complx  = "vector",
                        loglik  = "vector", 
                        hamm    = "vector",
                        hammexp = "vector",
                        tp = "vector",
                        fp = "vector",
                        fn = "vector",
                        pr = "vector",
                        sp = "vector",
                        sn = "vector",
                        fscore    = "vector",
                        skel.tp   = "vector",
                        skel.fp   = "vector",
                        skel.fn   = "vector",
                        order.fp  = "vector",
                        order.fn  = "vector", 
                        markov.fp = "vector",
                        markov.fn = "vector",
                        KLdist    = "vector",
                        time      = "numeric"
			),
         package = "sdnet"
	   )

setGeneric("cnFind", 
          function(object, complx=0, alpha=0,factor=1)
           standardGeneric("cnFind")
           )

setGeneric("cnParHist", 
          function(object)
           standardGeneric("cnParHist")
           )

setGeneric("cnProcTime", 
          function(object)
           standardGeneric("cnProcTime")
           )

#########################################################################
## for large sets of dags 

setClass("dagEvaluate",  
         representation(
			version = "character",
                        numnodes = "integer",
			nodes = "vector",
			numsamples = "integer", 
			maxcats="integer",
			cats="list",  
			maxpars="integer",  
                        parSlots="list",
			parLogliks="list",
			parComplx="list",
			parSampleSize="list",
                        numDags="integer", 
                        numPars="list",
                        loglik="vector",
			complx="vector",
			time = "numeric"), 
         package = "sdnet"
	   )

setGeneric("cnCatnetFromDagEvaluate", 
          function(object, index=1)
           standardGeneric("cnCatnetFromDagEvaluate")
           )

setGeneric("cnFindBnkl", 
          function(object, factor=NULL)
           standardGeneric("cnFindBnkl")
           )

#########################################################################
## Complete Partially Directed Acyclic Graphs

setClass("CPDAG",  
         representation(
                        numnodes = "integer",
                        nodes = "vector", 
                        edges = "list"
                        ),
         package = "sdnet"
	   )

setGeneric("dag2cpdag", 
           function(object)
           standardGeneric("dag2cpdag")
           )

#########################################################################
## Some constants

EPSDIFF <- exp(-16)

