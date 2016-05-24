print.gexf <- function(x, file=NA, replace=F, ...) {
################################################################################
# Printing method
################################################################################
  if (is.na(file)) {
    cat(x$graph)
  }
  else {
    output <- file(description=file,open="w",encoding='UTF-8')
    write(x$graph, file=output,...)
    close.connection(output)
    message('GEXF graph successfully written at:\n',normalizePath(file))
  }
}

summary.gexf <- function(object, ...) {
  ################################################################################
  # Printing method
  ################################################################################
  result <- list("N of nodes"=NROW(object$nodes), 
                 "N of edges"=NROW(object$edges),
                 "Node Attrs"=head(object$atts.definitions$node.att),
                 "Edge Attrs"=head(object$atts.definitions$edge.att))
  #class(result) <- "table"
  cat("GEXF graph object\n")
  result
}

.build.and.validate.gexf <- function(
  meta=list(creator="NodosChile", description="A graph file writing in R using \'rgexf\'",keywords="gexf graph, NodosChile, R, rgexf"),
  mode=list(defaultedgetype="undirected", mode="static"),
  atts.definitions=list(nodes = NULL, edges = NULL),
  nodesVizAtt=NULL,
  edgesVizAtt=NULL,
  nodes,
  edges,
  graph
  ) {
  
  # Shorcuts
  .c <- function(x,s) {
      if (NROW(x)) !all(colnames(x) %in% s)
      else return(FALSE)
  }
  .n <- function(x,s) !all(names(x) %in% s)
  
  # Check meta-data
  if (.n(meta,c("creator", "description","keywords")))
    stop("Invalid meta: Check names")
  else {
    for (i in meta) {
      if (!inherits(i,"character")) stop("Invalid meta: only character allowed")
    }
  }
  
  # Check mode
  
  if (.n(mode, c("defaultedgetype", "mode")))
    stop("Invalid mode: Check names")
  else {
    for (i in mode) {
      if (!inherits(i,"character")) stop("Invalid mode: only character allowed")
    }
  }
  
  # Checking atts definition
  if (.n(atts.definitions, c("nodes","edges"))) 
    stop("Invalid atts.definitions: Check names")
  else {
    for (i in atts.definitions) {
      
      # If its empty, then continue
      if (!length(i)) next
      
      if (!inherits(i,"data.frame")) stop("Invalid atts.definitions: only data-frames allowed")
      else {
        if (.c(atts.definitions$nodes, c("id","title","type")))
          stop("Invalid atts.definitions: Check -nodes- colnames")
        if (.c(atts.definitions$edges, c("id","title","type")))
          stop("Invalid atts.definitions: Check -nodes- colnames")
      }
    }
  }
  
  # Checking nodesVizAtt definition
  if (.n(nodesVizAtt, c("color","position","size","shape", "image")))
    stop("Invalid nodesVizAtt: Check names")
  else {    
    for (i in names(nodesVizAtt)) {
      if (i == "color" & .c(nodesVizAtt[[i]], c("r","g","b","a"))) 
        stop("Invalid nodesVizAtt: Check -color- colnames")
      else if (i == "position" & .c(nodesVizAtt[[i]], c("x","y","z")))
        stop("Invalid nodesVizAtt: Check -position- colnames")
      else if (i == "size" & .c(nodesVizAtt[[i]], "value"))
        stop("Invalid nodesVizAtt: Check -size- colname")        
      else if (i == "shape" & .c(nodesVizAtt[[i]], "value"))
        stop("Invalid nodesVizAtt: Check -shape- colname")
      else if (i == "image" & .c(nodesVizAtt[[i]], c("value","uri"))) 
        stop("Invalid nodesVizAtt: Check -image- colname")
    }
  }
  
  # Checking edgesVizAtt definition
  if (.n(edgesVizAtt, c("color","size","shape")))
    stop("Invalid edgesVizAtt: Check names")
  else {
    for (i in names(edgesVizAtt)) {
      if (i == "color" & .c(edgesVizAtt[["i"]], c("r","g","b","a"))) 
        stop("Invalid edgesVizAtt: Check -color- colnames")
      else if (i == "size" & .c(edgesVizAtt[["i"]], "value"))
        stop("Invalid edgesVizAtt: Check -size- colnames")
      else if (i == "shape" & .c(edgesVizAtt[["i"]], "value"))
        stop("Invalid edgesVizAtt: Check -shape- colname")
    }
  }
  
  # Checking nodes
  if (NROW(nodes)) {
    if (!all(c("id", "label") %in% colnames(nodes)))
      stop("Invalid nodes: Check colnames")
    else {
      nodes <- nodes[,unique(c("id","label",colnames(nodes)))]
      nodes <- nodes[,!grepl("^viz\\.[a-z]*\\.",colnames(nodes))]
    }
  }
  
  # Checking edges
  if (NROW(edges)) {
    if (!all(c("id", "source", "target", "weight") %in% colnames(edges)))
      stop("Invalid edges: Check colnames")
    else {
      edges <- edges[,unique(c("id", "source", "target", "weight",colnames(edges)))]
      edges <- edges[,!grepl("^viz\\.[a-z]*\\.",colnames(edges))]
    }
  }  
  
  # Returns the output
  output <- list(
    meta=unlist(meta),
    mode=unlist(mode),
    atts.definitions=atts.definitions,
    nodesVizAtt=nodesVizAtt,
    edgesVizAtt=edgesVizAtt,
    nodes=nodes,
    edges=edges,
    graph=graph
    )
  
  class(output) <- "gexf"
  
  return(output)
}

