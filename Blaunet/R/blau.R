blau <-
function(square.data, graph = NULL, directed.el = FALSE, node.ids = NULL, weights = NULL, ecology.ids = NULL, exclude = NULL, dimensions = NULL, memberships = NULL, complete.cases = FALSE){

  #exclude extraneous columns up front
  if (!is.null(exclude)){
    excluded <- correctFormat(exclude, square.data)
    square.data <- square.data[-excluded]    
  }

  #need: put this after everything has been created
  #right now can cause unexpected behavior
  if (complete.cases == TRUE){
    square.data <- square.data[complete.cases(as.data.frame(square.data)),]
  }
  else {
    square.data <- as.data.frame(square.data)
  }

  #now we need checks. the function checks for the non-nullity of each argument

  #blau object
  blauObj <- list() #should have a dataframe/list/matrix/etc for each option
  class(blauObj) <- 'blau'

  #ERROR CHECKS: it's vital that if the program is extended and new error checks for input format are needed that they be added here. The reason is simple: the input options are cleaned up and checks are displayed IMMEDIATELY to the user. There should be no waiting 60 seconds only to find an error in the input arguments. ALSO: getting errors out of the way and cleaning up the options arguments makes the following code MUCH easier to write and read.

  #checks whether arguments that should be length 1 are length 1
  if (!isCorrectLength(node.ids) || !isCorrectLength(ecology.ids) || !isCorrectLength(weights)) {print('Error in Argument Length')}

  #checks whether arguments are are in numeric form. if they're not, converts to numeric form. all column identifiers should be nuemric after this point.
  #if column names are needed (for writing, say), use colnames(OBJECT[colnumber])
  node.ids <- correctFormat(node.ids, square.data)
  ecology.ids <- correctFormat(ecology.ids, square.data)
  weights <- correctFormat(weights, square.data)
  dimensions <- correctFormat(dimensions, square.data)
  memberships <- correctFormat(memberships, square.data)

  #sorts the data frame by ecologies
  #sorting is IMPORTANT becuase we'd like all of the data to be grouped by ecologies.
  if (!is.null(ecology.ids)) {
    square.data <- square.data[order(square.data[, ecology.ids]), ]
  }

  
  #list of binary and continuous columns
  continuousCols <- NULL
  binaryCols <- NULL
  numericCols <- sapply(square.data, is.numeric)
  charCols <- which(numericCols == FALSE)


  for (colCyc in 1:ncol(square.data)) {    
    if (!isBinary(square.data[,colCyc]) && numericCols[colCyc] == TRUE) {
      continuousCols <- c(continuousCols, colCyc)
    }
    else if (isBinary(square.data[,colCyc]) && numericCols[colCyc] == TRUE) {
      binaryCols <- c(binaryCols, colCyc)
    }
  }


  #SINGLE ARGUMENT ASSIGNMENTS

  #idCol-- if there's no idCol, just use the row numbers. both names and numeric (and mixed) are good inputs
  if (is.null(node.ids)){ 
    tempNodeId <- c(1:nrow(square.data)) 
  }
  else { 
    tempNodeId <- as.vector(square.data[node.ids]) 
  }

  #ecologyId-- if no ecology.ids, everyone is in same ecology (#1). Else, people are placed in ecologies. 
  #this is like schoolID in the original program. 
  if (is.null(ecology.ids)) { 
    tempEcologyId <- rep(1, nrow(square.data))
  }
  else { 
    tempEcologyId <- as.vector(square.data[ecology.ids])
  }

  #put node and ecology identifiers together into one object
  blauObj$ids <- as.data.frame(cbind(tempNodeId, tempEcologyId))
  colnames(blauObj$ids) <- c('nodeId', 'ecologyId')

  #weights
  if (is.null(weights)) { 
    blauObj$weights <- as.matrix(rep(1, nrow(square.data)))
  } #default is a matrix of 1's
  else { 
    blauObj$weights <- as.matrix(square.data[weights])
  }


  #checks whether the graph argument is usable
  #if yes, puts it in a memory-efficient edgelist from the network package
  #the user MUST SPECIFY node ids in object that is turned into an edgelist
  #otherwise, there is no way we can assure that nodes are matched correctly
  if (!is.null(graph)) {
    if (class(graph) == 'network'){
      blauObj$graph <- graph
    }
    else {
      blauObj$graph <- network(graph, directed=directed.el)

      #make sure there are no nodes in the network that aren't in node ids
      for (name in network.vertex.names(blauObj$graph)){
        if (!any(blauObj$ids[,1] == name)){
          print(sprintf('Graph vertex with name %s is not present in node.ids.', name))
        }
      }
    }
  }


  #MULTIPLE ARGUMENT ASSIGNMENTS WITH DEFAULT SETTINGS

  #blauDimensions--subset the columns that are NOT already assigned AND are continuous variables 
  if (is.null(dimensions)) { 
    ignoredCols <- unique(c(node.ids, ecology.ids, weights, binaryCols, memberships, charCols))
    ignoredCols <- ignoredCols[!is.na(ignoredCols)]
    totalCols <- c(1:ncol(square.data))
    specifiedCols <- totalCols[-ignoredCols]

    blauObj$dimensions <- as.matrix(square.data[specifiedCols])
    } 
  else { #if not null, take specified columns and raise an error if there's overlap with columns reserved by other options
    ignoredCols <- unique(c(node.ids, ecology.ids, weights, binaryCols, memberships, charCols))
    if (length(intersect(dimensions, ignoredCols)) > 0) { 
      print('You have overlaps between specified Blau dimensions and other columns.')
    }
    else { 
      blauObj$dimensions <- as.matrix(square.data[dimensions])
    }
  } 

  #memberships-- just like with blauDimensions, if NULL, we automatically assign all binary unassigned variables to this category. if not NULL, we make sure there's no overlap and just take the user specified columns.
  if (is.null(memberships)){
    ignoredCols <- unique(c(node.ids, ecology.ids, weights, dimensions,continuousCols, charCols))
    ignoredCols <- ignoredCols[!is.na(ignoredCols)]
    totalCols <- c(1:ncol(square.data))
    specifiedCols <- totalCols[-ignoredCols]

    blauObj$memberships <- as.matrix(square.data[specifiedCols])
  }
  else { 
    ignoredCols <- unique(c(node.ids, ecology.ids, weights, dimensions,continuousCols, charCols))
    if (length(intersect(memberships,ignoredCols)) > 0) { 
      print('You have overlaps specified between membership columns and other columns.')
    }
    else { 
      blauObj$memberships <- as.matrix(square.data[memberships])
    }
  }


  #name the rows with the id names
  rownames(blauObj$ids) <- blauObj$ids[,1]
  rownames(blauObj$dimensions) <- blauObj$ids[,1]
  rownames(blauObj$memberships) <- blauObj$ids[,1]
  rownames(blauObj$weights) <- blauObj$ids[,1]

  if (!is.null(blauObj$primaryMembership)){
    rownames(blauObj$primaryMembership) <- blauObj$ids[,1]
  }

  #missing weight values
  presentObs <- complete.cases(blauObj$weights)

  blauObj$ids <- blauObj$ids[presentObs, , drop=FALSE]
  blauObj$dimensions <- blauObj$dimensions[presentObs, , drop=FALSE]
  blauObj$memberships <- blauObj$memberships[presentObs, , drop=FALSE]
  blauObj$weights <- blauObj$weights[presentObs, , drop=FALSE]

  if (!is.null(blauObj$primaryMembership)){
    blauObj$primaryMembership <- blauObj$primaryMembership[presentObs, , drop=FALSE]
  }

  #for the soul who decides to input a character matrix
  #this is here because datatypes in R can get confusing when both characters and numbers are stored in a data.frame
  #also, we're programming for potential R neophytes
  if (is.character(blauObj$dimensions)){
    print('The dimensions contain at least one character column. Dmensions must be numeric')
  }
  if (is.character(blauObj$memberships)){
    print('The memberships contain at least one character column. Memberships must be numeric')
  }
  if (is.character(blauObj$weights)){
    print('The weights contain at least one character element. Weights must be numeric')
  }

  #initilize null items for checks in subsequent functions
  #any new elements that are added to the blauobject should be intialized here
  blauObj$isInNiche <- NULL
  blauObj$topbounds <- NULL
  blauObj$lowbounds <- NULL
  blauObj$nodalLocal <- NULL
  blauObj$nodalGlobal <- NULL
  blauObj$nodalNetwork <- NULL
  blauObj$dyadic <- NULL

  #returns the data object
  return(blauObj)
}
