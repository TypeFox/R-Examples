require(RnavGraph) || stop("RnavGraph library not available")
local({
  ## Import the data
  ng.iris <- ng_data(name = "iris", data = iris[,1:4],
                     shortnames = c('s.L', 's.W', 'p.L', 'p.W'),
                     group = iris$Species,
                     labels = substr(iris$Species,1,2))
  
  ## get the variable graph node names
  V <- shortnames(ng.iris)
  
  ## create the linegraph and its complement
  G <- completegraph(V)
  LG <- linegraph(G)
  LGnot <- complement(LG)
  
  ## geberate NG_graph objects
  ng.lg <- ng_graph(name = '3D Transition', graph = LG, layout = 'circle')
  ng.lgnot <- ng_graph(name = '4D Transition', graph = LGnot, layout = 'circle')
  
  
  
  
  ## Create New Visualization Instructions
  ## Class
  setClass(
           Class="testVizClass",
           representation =
           representation(
                          ## No additional slots
                          ),
           contains = "NG_Visualization"
           )
  ## object creater
  myViz <- function(data,graph) {
    if(is(data,"NG_data") == FALSE){
      stop("data is no NG_data object.\n")
    }
    if(is(graph,"NG_graph") == FALSE){
      stop("graph is no NG_graph object.\n")
    }
    
    return(new(
               "testVizClass",
               graph = graph@name,
               data = data@name
               ))
  }
  
  ## methods
  setMethod(
            f = "initializeViz",
            signature = "testVizClass",
            definition = function(viz,ngEnv){
              print(paste('You switched to the graph', viz@graph))
              return(viz)
            })
  
  setMethod(
            f = "updateViz",
            signature = "testVizClass",
            definition = function(viz,ngEnv){
              print(paste('You current state is:', 
                          ngEnv$bulletState$from, 'to',
                          ngEnv$bulletState$to, 'and',
                          floor(ngEnv$bulletState$percentage*100), 'percent in between'))
              return(viz)
            })
  setMethod(
            f = "closeViz",
            signature = "testVizClass",
            definition = function(viz,ngEnv){
              print(paste('Bye Bye', viz@graph))
              return(viz)
            })
  

  ## custom visualization instructions
  vizNew <- myViz(ng.iris,ng.lg)
  
  ## visualization instructions for 2d scatterplots
  viz3dTransition <- ng_2d(ng.iris,ng.lg, glyphs=c("s.L","s.W","p.L","p.W"))
  viz4dTransition  <- ng_2d(ng.iris,ng.lgnot, glyphs=c("s.L","s.W","p.L","p.W"))
 
  
  ## pack them into list
  viz <- list(viz3dTransition, viz4dTransition, vizNew)	
  graphs <- list(ng.lg, ng.lgnot)
  
  ## start navGraph
  nav <- navGraph(data = ng.iris, graph = graphs, viz = viz, settings=list(tk2d=list(linked=FALSE)))
})

cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_iris.R", package="RnavGraph"),"\n\n\n"))
