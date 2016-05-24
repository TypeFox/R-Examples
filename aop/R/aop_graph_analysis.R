#convert_aop_to_graph function --------------------------------------
#' Convert AOP to Graph
#' 
#' 
#' Converts an AOP (encoded as \code{aop_cytoscape} object) to a 
#' \code{graphNEL} object.
#' 
#' This function converts an aop_cytoscape object to a \code{graphNEL} object.
#' This allows us to perform graph-based analyses of the AOP.
#' 
#' @param aop an object of class \code{aop_cytoscape}.
#' 
#' @return aop_graph a \code{graphNEL} object representation of the AOP
#' 
#' @examples
#' library(graph)
#' steatosis_json_file <- system.file("extdata", "steatosis_aop_json.cyjs",
#' package = "aop")
#' steatosis_aop <- convert_cytoscape_to_aop(steatosis_json_file)
#' steatosis_aop_graph <- convert_aop_to_graph(steatosis_aop)
#' plot(steatosis_aop_graph)
#' 
#' @import graph
#' @import Rgraphviz
#' @import methods
#' 
#' @export
convert_aop_to_graph <- function(aop){
  #Create the list of nodes
  #Then create a new graphNEL object
  node_graph_list <- names(aop@nodes)
  aop_graph <- new("graphNEL", nodes=node_graph_list, edgemode="directed")
  
  #Go through each edge and add it to the aop_graph object
  for(i in 1:length(aop@edges)){
    edge_source <- aop@edges[[i]]$data$source
    edge_target <- aop@edges[[i]]$data$target
    aop_graph <- addEdge(from=edge_source, to=edge_target, graph=aop_graph)
  }
  
  #return the aop as a graphNEL object
  return(aop_graph)
}

#---------------------------------------------------------------------

#Convert cytoscape to aop --------------------------------------------
#' Convert Cytoscape Graph to an AOP
#' 
#' Converts a cytoscape JSON file to an \code{aop_cytoscape-class} object.
#' 
#' This function converts a JSON file exported from Cytoscape into a
#' \code{aop_cytoscape-class} object. Once an \code{aop_cytoscape-class}
#' object, we can perform conversion to a \code{graphNEL} object, and then
#' perform graph-based analyses.
#' 
#' @param file a Cytoscape JSON file.
#' 
#' @return aop a \code{aop_cytoscape-class} object.
#' 
#' @examples
#' steatosis_json_file <- system.file("extdata", "steatosis_aop_json.cyjs",
#' package = "aop")
#' steatosis_aop <- convert_cytoscape_to_aop(steatosis_json_file)
#' 
#' @import rjson
#' 
#' @export
convert_cytoscape_to_aop <- function(file){
  aop_json <- fromJSON(file=file)
  
  #rename the list elements to be the $data$id
  node_ids <- NULL
  for(i in 1:length(aop_json$elements$nodes)){
    node_ids <- c(node_ids, aop_json$elements$nodes[[i]]$data$id)
  }
  names(aop_json$elements$nodes) <- node_ids
  node_list <- aop_json$elements$nodes
  
  edge_ids <- NULL
  for(i in 1:length(aop_json$elements$edges)){
    edge_ids <- c(edge_ids, aop_json$elements$edges[[i]]$data$id)
  }
  names(aop_json$elements$edges) <- edge_ids
  edge_list <- aop_json$elements$edges
  
  aop <- aop_cytoscape(name=aop_json$data$name, nodes=node_list, edges=edge_list)
  return(aop)
}


#AOP Backdoor Analysis --------------------------------------------
#' Backdoor Causal Network Analysis for AOPs
#' 
#' Performs a backdoor causal network analysis to idenitfy nodes/key events
#' which are sufficient to infer causality.
#' 
#' This function performs Pearl's backdoor analysis. Whereas Pearl was 
#' interested in identifying nodes which need to be measured to make a causal
#' statement, we are interested in identifying those nodes/key events which 
#' need to measured to say that an adverse outcome is likely to occur. It's 
#' essentially the same thing as Pearl, only a slightly different 
#' interpretation.
#' 
#' @param aop_graph a \code{graphNEL} object that encodes the AOP. Typically,
#' this would be output from the \code{convert_aop_to_graph function}.
#' 
#' @param ke_coord typically this is the molecular initiating event node, but
#' really, this is any node that you want as the starting/source point. For
#' instance, this is normally the point at which exposure to a stressor is
#' going to enter the AOP.
#' 
#' @param ao_coord typically this is the adverse outcome node.
#' 
#' @param measureable_nodes this param is not used yet. In the future this node
#' will be a \code{vector} of the nodes where an assay is available to measure
#' the node. In a future release this param will focus the backdoor algorithm
#' on finding only those nodes for which measurements can actually be taken, as
#' opposed to causal nodes regardless of our ability to measure them.
#' This allows for the assumption that AOP key events may or may not be 
#' measureable.
#'    
#' 
#' @return causal_nodes \code{vector} a vector of the names of the causal nodes. 
#' 
#' @examples
#' steatosis_json_file <- system.file("extdata", "steatosis_aop_json.cyjs",
#' package = "aop")
#' steatosis_aop <- convert_cytoscape_to_aop(steatosis_json_file)
#' steatosis_aop_graph <- convert_aop_to_graph(steatosis_aop)
#' aop_backdoor(steatosis_aop_graph, "391", "388")
#' 
#' @importFrom igraph graph.adjacency
#' @importFrom igraph get.all.shortest.paths
#' @importFrom igraph delete.vertices
#' @importFrom igraph get.adjacency
#' 
#' @export
aop_backdoor <- function(aop_graph, ke_coord, ao_coord, measureable_nodes=NULL){
  #Measurable nodes is a parameter that's not used yet, but will be in the future
  aop_adjm <- as(aop_graph, "matrix")
  
  #The ao_coord cannot be a child of the ke_coord
  if(as(aop_graph, "matrix")[ke_coord, ao_coord] == 1){
    return("ao_coord cannot be a direct descendent of the ke_coord")
  }
  
  #The ao_coord cannot be the same as the ke_coord
  if(ke_coord == ao_coord){
    return(print("The ke_coord and the ao_coord cannot be the same node"))
  }
  
  ##
  ##Not sure I want this next section to apply. I don't really mind if we have direct causal to be honest.
  #Get rid of all edges emanating from the ke_coord -- in Pearl this would be the exposure
  #Then make sure the graph and adjacency matrix agree
  #aop_adjm[which(rownames(aop_adjm) == ke_coord)] <- 0
  #aop_graph <- as(aop_adjm, "graphNEL")
  
  #Add edges between nodes that share a child and are measurable 
  for(i in 1:ncol(aop_adjm)){
    if(sum(aop_adjm[,i]) > 1){
      node_rows <- names(which(aop_adjm[,i] == 1))
      j <- 2
      while(j < length(node_rows) + 1){
        #Add an edge to the graph at the point where the nodes share the child node
        #I'm suppressing the warning that addEdge generates to inform you that it added the edge
        suppressWarnings(aop_graph <- addEdge(from=node_rows[j-1], to=node_rows[j], graph=aop_graph))
        j <- j+1
      }
    }
  }
  
  #Make sure the graph and adjacency matrix agree
  aop_adjm <- as(aop_graph, "matrix")  
  
  #Turn this into an igraph object and into an undirected graph
  aop_graph_u <- graph.adjacency(aop_adjm, mode="undirected") #graph
  aop_graph_m <- as(get.adjacency(aop_graph_u), "matrix")     #adjacency matrix
  
  #And now we do the iterative pathfinding -- it is important now to focus on identifying nodes and taking them out
  #We set a flag (continue_analysis) that will keep the while loop running as long as we get causal nodes showing up
  #We also create a vector, causal_nodes, to store the causal nodes as we find them
  continue_analysis <- TRUE
  causal_nodes <- NULL
  while(continue_analysis == TRUE){
    #this is the pathfinding algorithm from igraph -- it'll find the shortest path
    #if we remove the causal node after each pass, then we will find all of the backdoor paths
    causal_paths <- try(get.all.shortest.paths(aop_graph_u, from=ao_coord, to=ke_coord)$res)
    
    if(class(causal_paths) == "try-error"){
      return(print("The AO must occur after the KE"))
    }
    
    #as long as have some causal paths, the length will be greater than 0
    #we get the causal node, which is the second one in the list
    if(length(causal_paths) > 0){
      #causal_nodes <- c(causal_nodes, causal_paths[[1]][2])
      causal_nodes <- c(causal_nodes, rownames(get.adjacency(aop_graph_u))[causal_paths[[1]][2]])
      aop_graph_u <- delete.vertices(aop_graph_u, causal_paths[[1]][2])
    }
    else{
      continue_analysis <- FALSE
    }
  }
  #colnames(aop_adjm)[causal_nodes]
  #return(colnames(aop_adjm)[causal_nodes])
  return(causal_nodes)
}

