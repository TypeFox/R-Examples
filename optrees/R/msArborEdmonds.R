#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Arborescence Problems                                               #
#-----------------------------------------------------------------------------#

# msArborEdmonds --------------------------------------------------------------
#' Minimum cost arborescence with Edmonds' algorithm
#'
#' Given a connected weighted and directed graph, \code{msArborEdmonds} uses 
#' Edmonds' algorithm to find a minimum cost arborescence.
#' 
#' @details Edmonds' algorithm was developed by the mathematician and computer
#' scientist Jack R. Edmonds in 1967. Previously, it was proposed in 1965 by 
#' Yoeng-jin Chu and Tseng-hong Liu.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param source.node source node of the graph. It's node \eqn{1} by default.
#' @param stages.data logical value indicating if the function returns data of
#' each stage. Is \code{FALSE} by default.
#'
#' @return \code{msArborEdmonds} returns a list with:
#' tree.nodes vector containing the nodes of the minimum cost arborescence.
#' tree.arcs matrix containing the list of arcs of the minimum cost 
#' arborescence.
#' stages number of stages required.
#' stages.data complete data from each stage.
#' 
#' @references Chu, Y. J., and Liu, T. H., "On the Shortest Arborescence of a 
#' Directed Graph", Science Sinica, vol. 14, 1965, pp. 1396-1400.
#' 
#' Edmonds, J., "Optimum Branchings", Journal of Research of the National
#' Bureau of Standards, vol. 71B, No. 4, October-December 1967, pp. 233-240.
#' 
#' @seealso A more general function \link{getMinimumSpanningTree}.

msArborEdmonds <- function(nodes, arcs, source.node = 1, stages.data = FALSE) {
  
  # 0. Previous check
  # Check validity of the graph
  validGraph <- checkGraph(nodes, arcs, source.node, directed = TRUE)
  if (validGraph == FALSE) {
    # If graph is not valid stop algorithm before start
    stop("Invalid graph")
  }
  
  # 1. Search for a minimum cost arborescence
  # Initialize
  actualStage <- 1  # start counter
  stages <- c()  # vector to store stages
  arbor <- FALSE  # condition
  
  # Iterate until found an arborescence or the graph is not valid
  while (arbor == FALSE && validGraph == TRUE) {
    
    # Create a new name for the stage and save it
    stageName <- paste(c("stage", actualStage), collapse = "")
    stages <- c(stages, stageName)
    
    # Get minimum cost arcs reaching each node (except source)
    mcArcs <- getMinCostArcs(nodes, arcs)
    
    # Check arborescence with minimum cost arcs
    stage <- checkArbor(nodes, mcArcs)
    
    if (stage$arborescence == TRUE) {
      # If found arborescence notify and stop
      messArbor <- "Arborescence founded"
      # Build a list with the data generated so far
      stageData <- list(nodes = nodes, arcs = arcs, mcArcs = mcArcs,
                        arborescence = stage$arborescence,
                        arborArcs = stage$arbor.arcs)
      # Assign stageData to new object stageName
      assign(stageName, stageData)
      # Stop algorithm
      arbor <- TRUE
      
    } else {
      # Continue algorithm if no arborescence
      
      # Substract minimum costs to the arcs
      cheapArcs <- getCheapArcs(nodes, arcs)
      
      # Consider only arcs with zero cost
      zeroArcs <- getZeroArcs(nodes, cheapArcs)
      
      # Search cycle in graph only with zero cost arcs
      stageCycle <- searchZeroCycle(nodes, zeroArcs)
      
      if (stageCycle$cycle == FALSE) {
        # If no cycle graph not valid (this should never happen)
        print("There is no cycle. Invalid graph")
        validGraph <- FALSE
        return(validGraph)
        
      } else {
        # There is a cycle
        cycleNodes <- stageCycle$cycle.nodes  # cycle to compact
        
        # Compact nodes of the cycle with cheapArcs in new graph
        newGraph <- compactCycle(nodes, cheapArcs, cycleNodes)
        
        # X. Store data from each stage
        stageData <- list(nodes = nodes, arcs = arcs, mcArcs = mcArcs,
                          arborescence = stage$arborescence,
                          arborArcs = stage$arbor.arcs,
                          cheapArcs = cheapArcs, zeroArcs = zeroArcs,
                          cycleNodes = stageCycle$cycle.nodes,
                          cycleArcs = stageCycle$cycle.arcs,
                          newNodes = newGraph$nodes, newArcs = newGraph$arcs,
                          superNode = newGraph$super.node,
                          matches = newGraph$matches)
        # Assign stageData to new object stageName
        assign(stageName, stageData)
        
        # Update stage and repeat
        actualStage <- actualStage + 1     
        # New start
        nodes <- newGraph$nodes  # new nodes
        arcs <- newGraph$arcs  # new list of arcs
        
      }
      
    }
    #print(stages)
    #print(arbor)
    
  }
  
  # Save stages if needed
  if (stages.data == TRUE) {
    stageData <- list()
    for (i in 1: length(stages)) {
      stageData[[i]] <- eval(parse(text = stages[i]))
    }
  }
  
  # 2. Rebuild arborescence
  # Number of stages needed
  numStages <- length(stages)
  
  # Start with last stage with arborescence and previous without it
  lastStage <- eval(parse(text = stages[numStages]))  # make object
  numStages <- numStages - 1
  prevStage <- eval(parse(text = stages[numStages]))  # make object
  
  while (numStages >= 1) {
    # Check stages until reach stage1
    lastStage <- stepbackArbor(before = prevStage, after = lastStage)
    numStages <- numStages - 1; numStages
    prevStage <- eval(parse(text = stages[numStages]))
    #print(numStages)
  }
  
  # Resultado final
  TArcs <- matrix(lastStage$arborArcs, ncol = 3)
  
  # Number of stages
  nStages <- length(stages)
  
  if (stages.data == TRUE) {
    output <- list("tree.nodes" = nodes, "tree.arcs" = TArcs,
                   "stages" = nStages, "stages.data" = stageData) 
  } else {
    output <- list("tree.nodes" = nodes, "tree.arcs" = TArcs,
                   "stages" = nStages)
  }
  return(output)
  
}
#-----------------------------------------------------------------------------#