edge.list <- function(x) {
################################################################################
# Translate a edgelist to two objects list (nodes + edges)
################################################################################
  objClass <- class(x)
  nEdges <- NROW(x)
  nCols <- NCOL(x) == 2
  
  if (objClass %in% c("matrix", "data.frame")) {
    
    if (nCols) {
      # If it is not a factor
      if (!is.factor(x)) x <- factor(c(x[,1], x[,2]))
      edges <- matrix(unclass(x), byrow=FALSE, ncol=2)
      colnames(edges) <- c("source","target")
      nodes <- data.frame(id=1:nlevels(x), label=levels(x), stringsAsFactors=F)
      
      edgelist <- list(nodes=nodes, edges=edges)
      
      return(edgelist)
    }
    else stop("Insuficcient number of columns (", nCols,")")
  }
  else stop(objClass, 
                  " class not allowed, try with a \"matrix\" or a \"data.frame\"")
}

.defAtt <- function(x, parent) {
################################################################################
# Prints the nodes and edges att definition
################################################################################
  apply(x, MARGIN=1,
        function(x, PAR) {
          newXMLNode(name="attribute", parent=PAR, attrs=x)
        }, PAR=parent)
}

.addAtts <- function(tmpatt, attvec, tmpdoc=NULL) {
################################################################################
# Builds app proper XML attrs statement to be parsed by parseXMLAndAdd
################################################################################
  tmpatt <- data.frame(
    "for"=paste("att",attvec,sep=""), 
    value=unlist(tmpatt, recursive=FALSE), check.names=FALSE
    )

  for (i in attvec) 
    tmpdoc <- c(tmpdoc, .writeXMLLine("attvalue", tmpatt[i, ]) , sep="") 

  paste(c("<attvalues>", tmpdoc, "</attvalues>"), sep="", collapse="")
}

.writeXMLLine <- function(type, obj, finalizer=TRUE) {
################################################################################
# Builds as character whatever XML line is needed
################################################################################
  paste("<", type, " " ,
        paste(colnames(obj)[!is.na(obj)],obj[!is.na(obj)], sep="=\"", collapse="\" "),
        ifelse(finalizer, "\"/>","\">"), sep="")
}

.addNodesEdges <- function(dataset, PAR, type="node", doc) {
################################################################################
# Prints the nodes and edges
################################################################################  
  
  n <- NROW(dataset)
  vec <- 1:n
  xvars <- colnames(dataset)
  
  noattnames <- xvars[grep("(^att[0-9])|(^viz[.])", xvars, invert=T)]
  datasetnoatt <- dataset[, noattnames, drop=FALSE]
  
  # Parsing user-define attributes
  if (attributes <- length(grep("^att", xvars)) > 0) {
    attnames <- colnames(dataset)[grep("^att", xvars)]
    att <- dataset[,attnames, drop=FALSE]
    attvec <- 1:length(attnames)
  }
  
  # Parsing VIZ attributes
  if ((vizattributes <- length(grep("^viz[.]", xvars)) > 0)) {
    vizattnames <- colnames(dataset)[grep("^viz[.]", xvars)]
    
    # Color atts
    if ((vizcolors <- any(grepl("^viz[.]color",vizattnames)))) {
      vizcol.df <- dataset[,grep("^viz[.]color[.]", vizattnames, value=TRUE)]
      colnames(vizcol.df) <- gsub("^viz[.]color[.]", "", colnames(vizcol.df))
    }
    
    # Pos att
    if ((vizposition <- any(grepl("^viz[.]position",vizattnames)))) {
      vizpos.df <- dataset[,grep("^viz[.]position[.]", vizattnames, value=TRUE), drop=FALSE]
      colnames(vizpos.df) <- gsub("^viz[.]position[.]", "", colnames(vizpos.df))
    }
    
    # Size att
    if ((vizsize <- any(grepl("^viz[.]size",vizattnames)))) {
      vizsiz.df <- dataset[,grep("^viz[.]size[.]", vizattnames, value=TRUE), drop=FALSE]
      colnames(vizsiz.df) <- gsub("^viz[.]size[.]", "", colnames(vizsiz.df))
    }
    
    # Shape att
    if ((vizshape <- any(grepl("^viz[.]shape",vizattnames)))) {
      vizshp.df <- dataset[,grep("^viz[.]shape[.]", vizattnames, value=TRUE), drop=FALSE]
      colnames(vizshp.df) <- gsub("^viz[.]shape[.]", "", colnames(vizshp.df))
    }
    
    # Image att
    if ((vizimage <- any(grepl("^viz[.]image",vizattnames)))) {
      vizimg.df <- dataset[,grep("^viz[.]image[.]", vizattnames, value=TRUE), drop=FALSE]
      colnames(vizimg.df) <- c("value", "uri")
    }
    
    # Thickness att
    if ((viztness <- any(grepl("^viz[.]size",vizattnames)))) {
      vizthk.df <- dataset[,grep("^viz[.]size[.]", vizattnames, value=TRUE), drop=FALSE]
      colnames(vizthk.df) <- gsub("^viz[.]size[.]", "", colnames(vizthk.df))
    }
  }

  # Free memory
  rm(dataset)
  
  # Loop if there are not any attributes
  if (!attributes && !vizattributes) {
    for (i in vec) {
      parseXMLAndAdd(.writeXMLLine(type, datasetnoatt[i,]),parent=PAR)
    }
    return(NULL)
  }
  
  # Loop if only there are attributes
  if (attributes && !vizattributes) {
    
    for (i in vec) {      
      # Adding directly
      parseXMLAndAdd(
        paste(.writeXMLLine(type, datasetnoatt[i,], finalizer=FALSE), 
              .addAtts(att[i,], attvec), # Builds atts definition
              "</",type,">",sep=""),
        parent=PAR)
    }
    return(NULL)
  }
  
  # Loop if there are attributes and viz attributes
  for (i in vec) {
    # Node/Edge + Atts 
    if (attributes) {
      tempnode0 <- paste(
        .writeXMLLine(type, datasetnoatt[i,], finalizer=FALSE),
        .addAtts(att[i,], attvec), sep="")
    }
    else tempnode0 <- .writeXMLLine(type, datasetnoatt[i,], finalizer=FALSE)
    
    # Viz Att printing
    # Colors
    if (vizcolors) {
      tempnode0 <- paste(tempnode0, .writeXMLLine("color", vizcol.df[i,]),
                         sep="")
    }
    # Position
    if (vizposition) {
      tempnode0 <- paste(tempnode0, .writeXMLLine("position", vizpos.df[i,]),
                         sep="")
    }
    # Size
    if (vizsize) {
      tempnode0 <- paste(tempnode0, .writeXMLLine("size", vizsiz.df[i,1, drop=FALSE]),
                         sep="")
    }
    # Shape
    if (vizshape) {
      tempnode0 <- paste(tempnode0, .writeXMLLine("shape", vizshp.df[i,1,drop=FALSE]),
                         sep="")
    }
    # Image
    if (vizimage) {
      tempnode0 <- paste(tempnode0, .writeXMLLine("shape", vizimg.df[i,]),
                         sep="")
    }
    parseXMLAndAdd(sprintf("%s</%s>",tempnode0, type), parent=PAR)
  }
  return(NULL)
}

write.gexf <- function(
  ################################################################################  
  # Prints the gexf file
  ################################################################################
  nodes,
  edges,
  edgesLabel=NULL,
  edgesId=NULL,
  edgesAtt=NULL,
  edgesWeight=NULL,
  edgesVizAtt = list(color=NULL, size=NULL, shape=NULL),
  nodesAtt=NULL,
  nodesVizAtt = list(color=NULL, position=NULL, size=NULL, shape=NULL, image=NULL),
  nodeDynamic=NULL,
  edgeDynamic=NULL,
  digits=getOption("digits"),
  output = NA,
  tFormat="double",
  defaultedgetype = "undirected",
  meta = list(creator="NodosChile", description="A graph file writing in R using \"rgexf\"",keywords="gexf graph, NodosChile, R, rgexf"),
  keepFactors = FALSE,
  encoding="UTF-8"
) {
  
  ##############################################################################
  # CLASS CHECKS AND OTHERS CHECKS
  
  # Nodes
  if (is.data.frame(nodes) | is.matrix(nodes)) {
    if (NCOL(nodes) != 2) stop("-nodes- should have two columns not ", NCOL(nodes))
  }
  else stop("Invalid object type: -nodes- should be a two column data.frame or a matrix")
  
  # Edges
  if (is.data.frame(edges) | is.matrix(edges)) {
    if (NCOL(edges) != 2) stop("-edges- should have two columns not ", NCOL(edges))
  }
  else stop("Invalid object type: -edges- should be a two column data.frame or a matrix")
  
  # Edges Label
  .parseEdgesLabel(edgesLabel, edges)
  
  # Parsing Edges Id
  edgesId <- .parseEdgesId(edgesId, edges)
  
  # Parsing Edges Att
  nEdgesAtt <- .parseEdgesAtt(edgesAtt, edges)
  
  # Parsing edges Weight
  .parseEdgesWeight(edgesWeight, edges)
  
  # Parsing edges Viz Att
  nEdgesVizAtt <- .parseEdgesVizAtt(edgesVizAtt, edges)
  
  # Nodes Att
  nNodesAtt <- .parseNodesAtt(nodesAtt, nodes)
  
  # Parsing nodes Viz Atts
  nNodesVizAtt <- .parseNodesVizAtt(nodesVizAtt, nodes)
  
  # Checking the number of digits
  if (!is.integer(digits)) stop("Invalid number of digits ",digits,
                                ".\n Must be a number between 0 and 22")
  fmt <- sprintf("%%.%gg", digits)
  
  # Dynamics
  dynamic <- c(FALSE, FALSE)

  if (length(nodeDynamic) > 0) {
    if (is.data.frame(nodeDynamic) | is.matrix(nodeDynamic)) {
      if (NROW(nodeDynamic) == NROW(nodes)) dynamic[1] <- TRUE
      else stop("Insuficient number of rows: -nodeDynamic- (",NROW(nodeDynamic), " rows) should have the same number of rows than nodes there are (", NROW(nodes),")")
    } else stop("Invalid object type: -nodeDynamic- should be a two columns data.frame or a matrix")
  }
  
  if (length(edgeDynamic) > 0) {
    if (is.data.frame(edgeDynamic) | is.matrix(edgeDynamic)) {
      if (NROW(edgeDynamic) == NROW(edges)) dynamic[2] <- TRUE
      else stop("Insuficient number of rows: -edgeDynamic- (",NROW(edgeDynamic), " rows) should have the same number of rows than edges there are (", NROW(edges),")")
    } else stop("Invalid object type: -edgeDynamic- should be a two columns data.frame or a matrix")
  }
  
  ##############################################################################
  # Strings
  old.strAF <- getOption("stringsAsFactors")
  options(stringsAsFactors = FALSE)
  
  if (!any(dynamic)) mode <- "static" else mode <- "dynamic"
  
  # Starting xml
  xmlFile <- newXMLDoc(addFinalizer=TRUE)
  gexf <- newXMLNode(name="gexf", doc = xmlFile)
  
  # gexf att
  
  newXMLNamespace(node=gexf, namespace="http://www.gexf.net/1.2draft")
  newXMLNamespace(
    node=gexf, namespace="http://www.gexf.net/1.1draft/viz", prefix="viz")
  newXMLNamespace(
    node=gexf, namespace="http://www.w3.org/2001/XMLSchema-instance",
    prefix="xsi"
  ) 
  
  xmlAttrs(gexf) <- c( 
    "xsi:schemaLocation" = "http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd",
    version=1.2)
  
  # graph
  xmlMeta <- newXMLNode(name="meta", 
                        attrs=list(lastmodifieddate=as.character(Sys.Date())), 
                        parent=gexf)
  newXMLNode(name="creator", meta$creator, parent=xmlMeta)
  newXMLNode(name="description", meta$description, parent=xmlMeta)
  newXMLNode(name="keywords", meta$keywords, parent=xmlMeta)
  
  xmlGraph <- newXMLNode(name="graph", parent=gexf)
  if (mode == "dynamic") {
    
    # Fixing time factors
    if (keepFactors)
      for(i in 1:2) {
        if (dynamic[1]) nodeDynamic[,i] <- as.numeric(nodeDynamic[,i])
        if (dynamic[2]) edgeDynamic[,i] <- as.numeric(edgeDynamic[,i])
      }
    else
      for(i in 1:2) {
        if (dynamic[1]) nodeDynamic[,i] <- as.character(nodeDynamic[,i])
        if (dynamic[2]) edgeDynamic[,i] <- as.character(edgeDynamic[,i])
      }
    
    strTime <- c(unlist(nodeDynamic),unlist(edgeDynamic))

    endTime <- strTime
    
    # Checking start and ends
    strTime <- min(strTime, na.rm=TRUE)
    endTime <- max(endTime, na.rm=TRUE)
    
    xmlAttrs(xmlGraph) <- c(mode=mode, start=strTime, end=endTime,
                            timeformat=tFormat, defaultedgetype=defaultedgetype)
        
    # Replacing NAs
    if (dynamic[1]) {
      nodeDynamic[is.na(nodeDynamic[,1]),1] <- strTime
      nodeDynamic[is.na(nodeDynamic[,2]),2] <- endTime
    }
    if (dynamic[2]) {
      edgeDynamic[is.na(edgeDynamic[,1]),1] <- strTime
      edgeDynamic[is.na(edgeDynamic[,2]),2] <- endTime
    }
  } else {
    xmlAttrs(xmlGraph) <- c(mode=mode, defaultedgetype=defaultedgetype)
  }
  
  datatypes <- matrix(
    c(
      "string", "character",
      "integer", "integer",
      "float", "double",
      "boolean", "logical"
    ), byrow=TRUE, ncol =2)
  
  # nodes att definitions
  if (nNodesAtt > 0) {
    TIT <- colnames(nodesAtt)
    TYPE <- unlist(lapply(nodesAtt, typeof))
    CLASS <- unlist(lapply(nodesAtt, class))
    
    # Checks for factors (factor replacing is done later)
    if (keepFactors) TYPE[CLASS == "factor"] <- "integer"
    else TYPE[CLASS == "factor"] <- "string"
    
    nodesAttDf <- data.frame(
      id = paste("att",1:nNodesAtt,sep=""), 
      title = TIT, 
      type = TYPE
    )
    
    # Fixing datatype
    for (i in 1:NROW(datatypes)) {
      nodesAttDf$type <- gsub(datatypes[i,2], datatypes[i,1], nodesAttDf$type)
    }
    
    xmlAttNodes <- newXMLNode(name="attributes", parent=xmlGraph)
    xmlAttrs(xmlAttNodes) <- c(class="node", mode="static")
    .defAtt(nodesAttDf, parent=xmlAttNodes)
    
  } 
  else {
    nodesAttDf <- NULL
  }
  
  # edges att
  if (nEdgesAtt > 0) {
    TIT <- colnames(edgesAtt)
    TYPE <- unlist(lapply(edgesAtt, typeof))
    CLASS <- unlist(lapply(edgesAtt, class))
    
    # Checks for factors (factor replacing is done later)
    if (keepFactors) TYPE[CLASS == "factor"] <- "integer"
    else TYPE[CLASS == "factor"] <- "string"
    
    edgesAttDf <- data.frame(
      id = paste("att",1:nEdgesAtt,sep=""), 
      title = TIT, 
      type = TYPE
    )
    
    # Fixing datatype
    for (i in 1:NROW(datatypes)) {
      edgesAttDf$type <- gsub(datatypes[i,2], datatypes[i,1], edgesAttDf$type)
    }
    
    xmlAttEdges <- newXMLNode(name="attributes", parent=xmlGraph)
    xmlAttrs(xmlAttEdges) <- c(class="edge", mode="static")
    .defAtt(edgesAttDf, parent=xmlAttEdges)
  } 
  else {
    edgesAttDf <- NULL
  }
  
  # nodes vizatt
  ListNodesVizAtt <- NULL
  if (nNodesVizAtt > 0) {
    
    # Cohersing into data.frames
    nodesVizAtt <- lapply(nodesVizAtt, as.data.frame)
    
    for (i in names(nodesVizAtt)) {
      tmpAtt <- nodesVizAtt[[i]]
      
      if (i == "color") {
        colnames(tmpAtt) <- paste("viz.color", c("r","g","b","a"), sep=".")
      }
      else if (i == "position") {
        colnames(tmpAtt) <- paste("viz.position", c("x","y","z"), sep=".")
      }
      else if (i == "size") {
        colnames(tmpAtt) <- "viz.size.value"
        tmpAtt[,1] <- sprintf(fmt, tmpAtt[,1])
      }
      else if (i == "shape") {
        colnames(tmpAtt) <- "viz.shape.value"
      }
      else if (i == "image") {
        tmpAtt <- data.frame(x=rep("image",NROW(nodes)), viz.image.uri=tmpAtt)
        colnames(tmpAtt) <- c("viz.image.value","viz.image.uri")
      }
      
      if (length(ListNodesVizAtt) == 0) ListNodesVizAtt <- tmpAtt
      else ListNodesVizAtt <- data.frame(ListNodesVizAtt, tmpAtt)
      
      # Saving changes
      colnames(tmpAtt) <- gsub(sprintf("viz.%s.",i),"", colnames(tmpAtt))
      nodesVizAtt[[i]] <- tmpAtt
    }
  }
  
  # edges vizatt
  ListEdgesVizAtt <- NULL
  if (nEdgesVizAtt > 0) {
    
    # Cohersing into data.frames
    edgesVizAtt <- lapply(edgesVizAtt, as.data.frame)
    
    for (i in names(edgesVizAtt)) {
      tmpAtt <- edgesVizAtt[[i]]
      
      if (i == "color") {
        colnames(tmpAtt) <- paste("viz.color", c("r","g","b","a"), sep=".")
      }
      else if (i == "size") {
        colnames(tmpAtt) <- "viz.size.value"
        tmpAtt[,1] <- sprintf(fmt, tmpAtt[,1])
      }
      else if (i == "shape") {
        colnames(tmpAtt) <- "value"
      }
      
      if (length(ListEdgesVizAtt) == 0) ListEdgesVizAtt <- tmpAtt
      else ListEdgesVizAtt <- data.frame(ListEdgesVizAtt, tmpAtt)
      
      # Saving changes
      colnames(tmpAtt) <- gsub(sprintf("viz.%s.",i),"", colnames(tmpAtt))
      edgesVizAtt[[i]] <- tmpAtt
    }
  }
  
  ##############################################################################
  # The basic char matrix definition  for nodes
  
  if (dynamic[1]) nodeDynamic <- as.data.frame(nodeDynamic)
  
  if (nNodesAtt > 0) nodesAtt <- data.frame(nodesAtt)
  
  for (set in c(nodeDynamic, nodesAtt, ListNodesVizAtt)) {
    try(nodes <- data.frame(nodes, set), silent=TRUE)
  }
  
  # Naming the columns
  attNames <- nodesAttDf["id"]
  if (!is.null(nodeDynamic)) tmeNames <- c("start", "end") else tmeNames <- NULL
  
  colnames(nodes) <- unlist(c("id", "label", tmeNames, attNames, colnames(ListNodesVizAtt)))
  
  # Fixing factors
  if (keepFactors) {
    tofix <- which(lapply(nodes, class) %in% "factor")
    if (length(tofix)) {
      warning("Factor variables will be stored as -numeric-.",
              "\nIf you don't want this behavior, set -keepFactors- as -FALSE-.")
      nodes[,tofix] <- lapply(nodes[,tofix,drop=FALSE], as.numeric)
    }
  }
  else {
    tofix <- which(lapply(nodes, class) %in% "factor")
    if (length(tofix))
      nodes[,tofix] <- lapply(nodes[,tofix,drop=FALSE], as.character)
  }
  
  # NODES
  xmlNodes <- newXMLNode(name="nodes", parent=xmlGraph)
  .addNodesEdges(nodes, xmlNodes, "node")
  
  ##############################################################################
  # The basic dataframe definition  for edges  
  
  if (dynamic[2]) edgeDynamic <- as.data.frame(edgeDynamic)
  
  if (nEdgesAtt > 0) edgesAtt <- data.frame(edgesAtt)
  
  # Adding edge id
  try(edgesId <- cbind(edgesId, edgesLabel), silent=TRUE)
  edges <- cbind(edgesId, edges)
  for (set in c(edgeDynamic, edgesAtt, ListEdgesVizAtt)) {
    try(edges <- data.frame(edges, set), silent=TRUE)
  }
  
  # Naming the columns
  attNames <- edgesAttDf["id"]
  if (!is.null(edgeDynamic)) tmeNames <- c("start", "end") else tmeNames <- NULL
  
  # Generating weights
  if (!length(edgesWeight))  edgesWeight <- 1
  edges <- data.frame(edges, x=as.numeric(edgesWeight))
  edges$x <- sprintf(fmt, edges$x)
  
  # Seting colnames
  if (length(edgesLabel) > 0) edgesLabelCName <- "label"
  else edgesLabelCName <- NULL
  colnames(edges) <- unlist(c("id", edgesLabelCName, "source", "target", 
                              tmeNames, attNames, colnames(ListEdgesVizAtt),
                              "weight"))
  
  # EDGES
  xmlEdges <- newXMLNode(name="edges", parent=xmlGraph)
  
  # Fixing factors
  if (keepFactors) {
    for (i in colnames(edges)) {
      if (class(edges[[i]]) == "factor") edges[[i]] <- as.numeric(edges[[i]])
    }
  }
  else {
    for (i in colnames(edges)) {
      if (class(edges[[i]]) == "factor") edges[[i]] <- as.character(edges[[i]])
    } 
  }
  
  # Adding edges
  .addNodesEdges(edges, xmlEdges, "edge")
  
  # Edges Label (for data frame)
  if (length(edgesLabel) == 0) edgesLabel <- edges[,"id"]
    
  results <- .build.and.validate.gexf(
    meta=meta,
    mode=list(defaultedgetype=defaultedgetype, mode=mode),
    atts.definitions=list(nodes = nodesAttDf, edges = edgesAttDf),
    nodesVizAtt=nodesVizAtt,
    edgesVizAtt=edgesVizAtt,
    nodes=nodes,
    edges=cbind(edges,edgesLabel),
    graph=saveXML(xmlFile, encoding=encoding)
    )
  
  # Strings As Factors
  options(stringsAsFactors = old.strAF)
  
  # Fixing 
  for (viz in c("color", "size", "shape", "position")) 
    results$graph <- gsub(sprintf("<%s",viz), sprintf("<viz:%s", viz), 
                          results$graph, fixed=TRUE)
  
  
  # Returns
  if (is.na(output)) {
    return(results)
  } else {
    print(results, file=output, replace=TRUE)
  }
}
