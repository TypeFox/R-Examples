################################################################################
# COLLECTION OF FUNCTIONS FOR EDGES ANALYSIS
################################################################################

check.dpl.edges <- function(edges, undirected=FALSE, order.edgelist=TRUE) {
################################################################################
# Checks for duplicated edges, and switchs between source and target
# (optionally).
################################################################################  
  srce <- edges[,1]
  trgt <- edges[,2]

  if (any(!is.finite(edges) | is.null(edges)))
    stop("No NA, NULL or NaN elements can be passed to this function.")
 
  nedges <- length(srce)
  
  result <- .C("RCheckDplEdges", 
     as.double(srce),           # Input Source
     as.double(trgt),           # Input Target
     as.integer(undirected),    # Tells the function if the graph is undirected
     "source" = as.double(      # Output Source
       vector("double", nedges) 
       ),
     "target" = as.double(      # Output Target
       vector("double", nedges)
     ),
     "reps" = as.double(        # Output Target
       vector("double", nedges)
     ), PACKAGE="rgexf"
     )
  
  result <- data.frame(source=result$source, target=result$target, 
                       reps=result$reps, check.names=FALSE)

  if (order.edgelist) 
    result <- result[order(result[,1], result[,2]),]
  
  return(result)
}

switch.edges <- function(edges) {
################################################################################
# Orders pairs of edges by putting the lowest id first as source
################################################################################
  if (any(is.na(edges) | is.null(edges) | is.nan(edges))) 
    stop("No NA, NULL or NaN elements can be passed to this function.")

  result <- .C(
    "RSwitchEdges",
    as.integer(NROW(edges)),
    as.double(edges[,1]),
    as.double(edges[,2]),
    "source" = as.double(              # Output Source
      vector("double", NROW(edges)) 
    ),
    "target" = as.double(              # Output Target
      vector("double", NROW(edges))
    ), PACKAGE="rgexf"
    )
  
  return(
    data.frame(
      source=result$source, 
      target=result$target,
      check.names=FALSE)
  )
}

#try(dyn.unload("src/RCheckDplEdges.so"))
#dyn.load("src/RCheckDplEdges.so")

#relations <- cbind(c(1,1,3,4,2,5,6), c(2,3,1,2,4,1,1))


#check.dpl.edges(relations)
#ordAndCheckDplEdges(relations, undirected=FALSE)
#switch.edges(relations)
