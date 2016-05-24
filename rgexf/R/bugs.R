checkTimes <- function(x, format='date') {
  ################################################################################
  # A a function of format, checks that all data has the correct format
  ################################################################################
  
  # Defining regular expressions to match
  if (format=='date') {
    match <- '^[0-9]{4}[-]{1}[01-12]{2}[-]{1}[01-31]{2}$'
  }
  else if (format == 'dateTime') {
    match <- '^[0-9]{4}[-][0-9]{2}[-][0-9]{2}[T][0-9]{2}[:][0-9]{2}[:][0-9]{2}$'
  }
  else if (format == 'float') {
    match <- '^[0-9]+[.]{1}[0-9]+$'
  }
  
  # Defining matchin function
  FUN <- function(x, pattern,...) {
    x <- grepl(pattern, x)
  }
  # listapply
  result <- lapply(x, FUN, pattern=match)
  result <- unlist(result, use.names=F)
  return(result)
}

#test <- c('2012-01-17T03:46:41', '2012-01-17T03:46:410')

#checkTimes(test, format='dateTime')
#checkTimes('2012-02-01T00:00:00', 'dateTime')

.parseNodesVizAtt <- function(nodesVizAtt, nodes) {
################################################################################
# Parses nodes viz attributes checking dimentions, classes and names
################################################################################
  if (any(lapply(nodesVizAtt, length) > 0)) {
    supportedNodesVizAtt <- c("color", "position", "size", "shape", "image")
    if (all(names(nodesVizAtt) %in% supportedNodesVizAtt)) {
      if (all(lapply(nodesVizAtt, NROW) == NROW(nodes))) {
        return(length(nodesVizAtt))
      }
      else {
        nodesVizAtt <- lapply(nodesVizAtt, NROW)
        nodesVizAtt <- nodesVizAtt[nodesVizAtt != NROW(nodes)]
        stop("Insuficient number of \"nodeVizAtt\" rows: The atts ",
             paste(names(nodesVizAtt), unlist(nodesVizAtt), sep=" (", collapse=" rows), "),")\n",
             "Every att should have the same number of rows than nodes there are (",NROW(nodes),")")
      }
    }
    else {
      noviz <- names(nodesVizAtt)
      noviz <- noviz[!(noviz %in% supportedNodesVizAtt)]
      stop("Invalid \"nodeVizAtt\": ",noviz,"\nOnly \"color\", \"position\", \"size\", \"shape\" and \"image\" are supported")
    }
  }
  else return(0)
}

.parseEdgesVizAtt <- function(edgesVizAtt, edges) {
################################################################################
# Parses edges viz attributes checking dimentions, classes and names
################################################################################
  if (any(lapply(edgesVizAtt, length) > 0)) {
    supportedEdgeVizAtt <- c("color", "size", "shape")
    if (all(names(edgesVizAtt) %in% supportedEdgeVizAtt)) {
      if (all(lapply(edgesVizAtt, NROW) == NROW(edges))) {
        return(length(edgesVizAtt))
      }
      else {
        edgesVizAtt <- lapply(edgesVizAtt, NROW)
        edgesVizAtt <- edgesVizAtt[edgesVizAtt != NROW(edges)]
        stop("Insuficient number of \"edgeVizAtt\" rows: The atts ",
             paste(names(edgesVizAtt), unlist(edgesVizAtt), sep=" (", collapse=" rows), "),")\n",
             "Every att should have the same number of rows than edges there are (",NROW(edges),")")
      }
    }
    else {
      noviz <- names(edgesVizAtt)
      noviz <- noviz[!(noviz %in% supportedEdgeVizAtt)]
      stop("Invalid \"edgesVizAtt\": ",noviz,"\nOnly \"color\", \"size\" and \"shape\" are supported")
    }
  }
  else return(0)
}

.parseEdgesWeight <- function(edgesWeight, edges) {
################################################################################
# Parses edges weights checking dimentions and classes
################################################################################
  if (length(edgesWeight) > 0) {
    if (is.vector(edgesWeight) | is.data.frame(edgesWeight) | is.matrix(edgesWeight)) {
      if (NROW(edgesWeight) != NROW(edges)) stop("\"edgesWeight\" should have the same number of rows than edges there are (", NROW(edges),")")
      if (NCOL(edgesWeight) > 1) stop("\"edgesWeight should have only one column\"")
    }
    else stop("Invalid object type: \"edgesWeight\" should be a one column data.frame, a matrix or a vector")
  }
}

.parseEdgesAtt <- function(edgesAtt, edges) {
################################################################################
# Parses edges attributes checking dimentions and classes
################################################################################
  if ((nEdgesAtt <- length(edgesAtt)) > 0) {
    if (is.data.frame(edgesAtt) | is.matrix(edgesAtt) | is.vector(edgesAtt)) {
      if (NROW(edgesAtt) != NROW(edges)) stop(paste("\"edgesAtt\" should have the same number of rows than edges there are (", NROW(edges),")",sep=""))
      else return(nEdgesAtt)
    }
    else stop("Invalid object type: \"edgesAtt\" should be a data.frame, a matrix or a vector")
  }
  else return(0)
}

.parseEdgesId <- function(edgesId, edges) {
################################################################################
# Parses edges Ids and if does not exists it creates them
################################################################################
  if (length(edgesId) > 0) {
    if (is.data.frame(edgesId) | is.matrix(edgesId) | is.vector(edgesId)) {
      if (NCOL(edgesId) != 1) stop("\"edgesId\" should have one column not ", NCOL(edgesId))
    }
    else stop("Invalid object type: \"edgesId\" should be a one column data.frame or a matrix")
  }
  else return(data.frame(id=0:(NROW(edges) - 1)))
}

.parseNodesAtt <- function(nodesAtt, nodes) {
################################################################################
# Parses nodes attributes checking dimentions
################################################################################
  if ((nNodesAtt <- length(nodesAtt)) > 0) {
    if (is.data.frame(nodesAtt) | is.matrix(nodesAtt) | is.vector(nodesAtt)) {
      if (NROW(nodesAtt) != NROW(nodes)) stop("Insuficient number of rows: \"nodesAtt\" (", NROW(nodesAtt)," rows) should have the same number of rows than nodes there are (", NROW(nodes),")")
      else return(nNodesAtt)
    }
    else stop("Invalid object type: \"nodesAtt\" should be a data.frame, a matrix or a vector")
  }
  else return(0)
}

.parseEdgesLabel <- function(edgesLabel, edges) {
################################################################################
# Parses edges labels checking dimentions
################################################################################
  if (length(edgesLabel) > 0) {
    if (is.data.frame(edgesLabel) | is.matrix(edgesLabel) | is.vector(edgesLabel)) {
      if (NCOL(edgesLabel) != 1) stop("\"edgesLabel\" should have one column not ", NCOL(edgesLabel))
    }
    else stop("Invalid object type: \"edgesLabel\" should be a one column data.frame or a matrix")
  }
}

.parseDataTypes <- function(x, keepFactors=TRUE) {
################################################################################
# Parses edges labels checking dimentions
################################################################################
  # Whether to keep factors as numeric values or not
  if (keepFactors) x <- as.numeric(x)
  else x <- as.character(x)
  
  # Data
  type <- typeof(x)
  if (type == "character") return("string")
  else if (type == "double") return("float")
  else if (type == "logical") return("boolean")
  else return(type)        
}
