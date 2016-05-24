#  File networkDynamic/R/import.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#scripts for loading in dynamic network file formats



#load a .son format file
#many features not supported yest
read.son <- function(file,guess.TEA=TRUE){
  
  alphaIDs<-NULL
  #find the index of the header line
  nodeStartLine <- grep("^NodeId",readLines(file))
  if (length(nodeStartLine) < 1){
    #try the other alternate
    nodeStartLine <- grep("^AlphaId",readLines(file))
    alphaIDs<-TRUE
  }
  if (length(nodeStartLine)<1){
    stop("Unable to locate the header line beginning with 'NodeId' or 'AlphaId'. Not a correctly formatted .son file")
  }
  #find the index of the edges split
  hasArcs<-TRUE
  arcStartLine <- grep("FromId",readLines(file))
  fileLength<-length(readLines(file))
  if (length(arcStartLine)<1){
    warning("Unable to locate an arc header line containing 'FromId'. Either this .son file is not correctly formatted or the network contains no edges")
    # set arcStartLine to length of file
    arcStartLine<-fileLength+1
    hasArcs<-FALSE
  }
  # case: arcs header exists, but there are no arcs
  if (arcStartLine==fileLength){
    hasArcs<-FALSE
  }
  # find the index of the (possible clusters split)
  # clusters are not implemented, but we don't want the parser to break if they are there
  clusterStartLine<-grep("^ClusterId",readLines(file))
  nrows<- -1
  if (length(clusterStartLine)>0){
    warning('parsing of clusters is not currently implemented in read.son, clusters skipped')
    nrows<-clusterStartLine-(arcStartLine+1)
  }
  
  noderows <- read.table(file,header=TRUE,sep="\t",as.is=TRUE,skip=nodeStartLine-1,nrows=arcStartLine-(nodeStartLine+1),stringsAsFactors=FALSE)
  # only read arcs if there are arcs to read
  if (hasArcs){
    arcrows <- read.table(file,header=TRUE,sep="\t",as.is=TRUE,skip=arcStartLine-1,nrows=nrows,stringsAsFactors=FALSE)
  } else {
    # no arcs, so need a dummy data frame
    arcrows<-as.data.frame(matrix(numeric(0),nrow=0,ncol=4,dimnames=list(c(),c('FromId','ToId','StartTime','EndTime'))))
  }
  
  # if we are doing alphaIds, collect them
  if (!is.null(alphaIDs)){
    alphaIDs<-unique(noderows$'AlphaId')
  } else { # if not, verify numeric ids
    idSet<-unique(noderows$'NodeId')
    if(!all(idSet==1:length(idSet))){
      stop("'NodeId' values must be a set of integers from 1 to the size of the network")
    }
  }
  
  
  #figure out the order of the node column headings
  idIndex <-1
  nodeStartIndex <- which(names(noderows)=="StartTime")
  nodeEndIndex <- which(names(noderows)=="EndTime")
  if (length(nodeStartIndex)<1){
    stop("Unable to locate a node column for 'StartTime' to provide vertex onset times")
  }
  if (length(nodeEndIndex)<1){
    if (length(nodeStartIndex)>0){
      # use start time as end time
      message("Unable to locate a node column for 'EndTime' to provide vertex terminus times, using values from 'StartTime'")
      nodeEndIndex<-nodeStartIndex
    } else {
      stop("Unable to locate a node column for 'EndTime' to provide vertex terminus times")
    }
    
  }

  #figure out the order of arc column headings
  fromIdIndex <-1
  toIdIndex <- 2
  arcStartIndex <- which(names(arcrows)=="StartTime")
  arcEndIndex <- which(names(arcrows)=="EndTime")
  if (length(arcStartIndex)<1){
    warning("Unable to locate an arc column for 'StartTime' to provide arc/edge onset times")
  }
  if (length(arcEndIndex)<1){
    if(length(arcStartIndex)>0){
    warning("Unable to locate an arc column for 'EndTime' to provide arc/edge terminus times, using values from 'StartTime")
     arcEndIndex<-arcStartIndex
    } else {
      warning("Unable to locate an arc column for 'EndTime' to provide arc/edge terminus times")
    }
  }
  
  # if doing alpha ids, replace the alphas with numeric values for vertices and edge records
  if (!is.null(alphaIDs)){
    noderows$'AlphaId'<-match(noderows$'AlphaId',alphaIDs)
    arcrows[,fromIdIndex]<-match(arcrows[,fromIdIndex],alphaIDs)
    arcrows[,toIdIndex]<-match(arcrows[,toIdIndex],alphaIDs)
  }
  dnet<-networkDynamic(edge.spells=arcrows[,c(arcStartIndex,arcEndIndex,fromIdIndex,toIdIndex)],vertex.spells=noderows[,c(nodeStartIndex,nodeEndIndex,idIndex)])
  
  # process the vertex attributes
  # TODO: this is going to be slow, should do somehow inside the construction stage
  vertAttrs<-colnames(noderows)[!colnames(noderows)%in%c('AlphaId','NodeId','StartTime','EndTime')]
  for (vAttr in vertAttrs){
    # determine if the attribute values change
    attrIndex<-match(vAttr,colnames(noderows))
    if(guess.TEA && nrow(unique(noderows[,c(idIndex,attrIndex)]))==length(unique(noderows[,idIndex]))){
    # if they don't change, set as regular attribute
      set.vertex.attribute(dnet,attrname=vAttr,value=noderows[,attrIndex],v=noderows[,idIndex])
    # if the attribute is named 'Label' and not TEA, also set as vertex.names
      if (vAttr=='Label'){
        network.vertex.names(dnet)<-get.vertex.attribute(dnet,'Label')
      }
      
    } else {
    # if they do change, loop and set as TEA
    # this will be horribly slow  
      for (r in seq.int(nrow(noderows))){
        activate.vertex.attribute(dnet,prefix=vAttr,value=noderows[r,attrIndex],v=noderows[r,idIndex],onset=noderows[r,nodeStartIndex],terminus=noderows[r,nodeEndIndex])
      }
    }
  }
  
  # process the edge attributes
  edgeAttrs<-colnames(arcrows)[!colnames(arcrows)%in%c('FromId','ToId','StartTime','EndTime')]
  
  for (eAttr in edgeAttrs){
    # determine if the attribute values change
    attrIndex<-match(eAttr,colnames(arcrows))
    collapsedDyads<-unique(arcrows[,c(fromIdIndex,toIdIndex,attrIndex)])
    if(guess.TEA && nrow(collapsedDyads)==nrow(unique(arcrows[,c(fromIdIndex,toIdIndex)]))){
      # if they don't change, set as regular attribute
      # but still have to loop over edges to get eids :-(
      for (r in seq.int(nrow(collapsedDyads))){
        eid<-get.edgeIDs(dnet,v=collapsedDyads[r,1],alter=collapsedDyads[r,2])
        set.edge.attribute(dnet,attrname=eAttr,value=collapsedDyads[r,3],e=eid)
      }
    } else {
      # if they do change, loop and set as TEA
      # this will be horribly slow  
      for (r in seq.int(nrow(arcrows))){
        eid<-get.edgeIDs(dnet,v=arcrows[r,fromIdIndex],alter=arcrows[r,toIdIndex])
        activate.edge.attribute(dnet,prefix=eAttr,value=arcrows[r,attrIndex],e=eid,onset=arcrows[r,arcStartIndex],terminus=arcrows[r,arcEndIndex])
      }
    }
  }
  
  
  
  # if doing alphaIds, set vertex names and pid.
  if (!is.null(alphaIDs)){
   #if label is not specified, use id as label
   network.vertex.names(dnet)<-alphaIDs
   dnet%n%'vertex.pid'<-'vertex.names'
  }  
  
  return(dnet)
}
