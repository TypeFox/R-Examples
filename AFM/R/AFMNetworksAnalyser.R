require(igraph)

HASHSIZE<-512*512


setOldClass("igraph")

#' AFM image networks analysis class
#' 
#' A S4 class to handle the networks calculation 
#'
#' @slot originalGraph a list of \code{\link{igraph}}
#' @slot skeletonGraph a list of \code{\link{igraph}}
#' @slot heightNetworksslider used multiplier of heights to facilitate analysis
#' @slot filterNetworkssliderMin used filter minimum value to facilitate analysis
#' @slot filterNetworkssliderMax used filter maximum value to facilitate analysis
#' @slot libVersion version of the AFM library used to perform the analysis
#' @slot updateProgress a function to update a graphical user interface
#' @name AFMImageNetworksAnalysis-class
#' @rdname AFMImageNetworksAnalysis-class
#' @exportClass AFMImageNetworksAnalysis
#' @author M.Beauvais
AFMImageNetworksAnalysis<-setClass("AFMImageNetworksAnalysis",
                                   slots = c(originalGraph="igraph", 
                                             skeletonGraph="igraph",
                                             heightNetworksslider="numeric",
                                             filterNetworkssliderMin="numeric",
                                             filterNetworkssliderMax="numeric",
                                             libVersion="character",
                                             updateProgress="function"))

#' Constructor method of AFMImageNetworksAnalysis Class.
#'
#' @param .Object an AFMImageNetworksAnalysis Class
#' @param originalGraph a list of \code{\link{igraph}}
#' @param skeletonGraph a list of \code{\link{igraph}}
#' @param heightNetworksslider used multiplier of heights to facilitate analysis
#' @param filterNetworkssliderMin used filter minimum value to facilitate analysis
#' @param filterNetworkssliderMax used filter maximum value to facilitate analysis
#' @param libVersion version of the AFM library used to perform the analysis
#' @rdname AFMImageNetworksAnalysis-class
#' @export
setMethod("initialize", "AFMImageNetworksAnalysis", function(.Object, originalGraph, 
                                                             skeletonGraph,
                                                             heightNetworksslider,
                                                             filterNetworkssliderMin,
                                                             filterNetworkssliderMax,
                                                             libVersion)  
{
  if(!missing(originalGraph)) .Object@originalGraph<-originalGraph
  if(!missing(skeletonGraph)) .Object@skeletonGraph<-skeletonGraph
  if(!missing(heightNetworksslider)) .Object@heightNetworksslider<-heightNetworksslider
  if(!missing(filterNetworkssliderMin)) .Object@filterNetworkssliderMin<-filterNetworkssliderMin
  if(!missing(filterNetworkssliderMax)) .Object@filterNetworkssliderMax<-filterNetworkssliderMax
  if(!missing(libVersion)) .Object@libVersion<-libVersion
  #validObject(.Object)      
  return(.Object)
})


#' Wrapper function AFMImageNetworksAnalysis
#'
#' @rdname AFMImageNetworksAnalysis-class
#' @export
AFMImageNetworksAnalysis <- function() {
  return(new("AFMImageNetworksAnalysis"))
}

#' Calculate networks on the surface
#'
#' \code{calculateNetworks} update  \code{\link{AFMImageNetworksAnalysis}}
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param AFMImageNetworksAnalysis n \code{\link{AFMImageNetworksAnalysis}} to store the results of networks analysis
#' 
#' @name calculateNetworks
#' @rdname calculateNetworks-methods
#' @exportMethod calculateNetworks
#' @author M.Beauvais
setGeneric(name= "calculateNetworks", 
           def= function(AFMImageNetworksAnalysis, AFMImage) {
             return(standardGeneric("calculateNetworks"))
           })

#' @rdname calculateNetworks-methods
#' @aliases calculateNetworks,AFMImage-method
setMethod(f="calculateNetworks", "AFMImageNetworksAnalysis",
          definition= function(AFMImageNetworksAnalysis, AFMImage) {
            
            counter<-0
            totalLength<-2
            if (!is.null(AFMImageNetworksAnalysis@updateProgress)&&
                is.function(AFMImageNetworksAnalysis@updateProgress)&&
                !is.null(AFMImageNetworksAnalysis@updateProgress())) {
              text <- paste0("Creating networks")
              AFMImageNetworksAnalysis@updateProgress(value= 0, detail = text)
              
              counter<-counter+1
              value<-counter / totalLength
              text <- paste0("Creating networks", round(counter, 2),"/",totalLength)
              AFMImageNetworksAnalysis@updateProgress(value= value, detail = text)
              print("update")
            }
            
            AFMImageNetworksAnalysis@originalGraph<-calculateIgraph(AFMImageNetworksAnalysis= AFMImageNetworksAnalysis, AFMImage = AFMImage)
            
            if (!is.null(AFMImageNetworksAnalysis@updateProgress)&&
                is.function(AFMImageNetworksAnalysis@updateProgress)&&
                !is.null(AFMImageNetworksAnalysis@updateProgress())) {
              text <- paste0("Creating networks skeleton")
              AFMImageNetworksAnalysis@updateProgress(value= 0, detail = text)
              
              counter<-counter+1
              value<-counter / totalLength
              text <- paste0("Creating networks", round(counter, 2),"/",totalLength)
              AFMImageNetworksAnalysis@updateProgress(value= value, detail = text)
              print("update")
            }
            
            AFMImageNetworksAnalysis<-calculateNetworkSkeleton(AFMImageNetworksAnalysis= AFMImageNetworksAnalysis, AFMImage = AFMImage)
            
            return(AFMImageNetworksAnalysis)
          })


#' Get vertex id from x,y coordinates
#'
#' \code{getVertexId} return the vertexId
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param x coordinates in x axis
#' @param y coordinates in y axis
#' @author M.Beauvais
#' @export
getVertexId<-function(AFMImage,x,y) {
  if ((x<1)||(x>AFMImage@samplesperline)||
      (y<1)||(y>AFMImage@lines)) return(-1)
  #print(paste("getVertexId",x,y,as.numeric(x+HASHSIZE*y)))
  #return(as.numeric(x+AFMImage@samplesperline*y))
  return(as.numeric(x+HASHSIZE*y))
  
}

#' Get x,y coordinates from vertex id
#'
#' \code{getCoordinatesFromVertexId} return a list x,y coordinates
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param vId the vertex id
#' @author M.Beauvais
#' @export
getCoordinatesFromVertexId<-function(AFMImage, vId) {
  vertexId<-as.numeric(vId)
  y<-floor(vertexId/HASHSIZE)
  x<-vertexId-y*HASHSIZE
  return(c(x,y))
}

#' Get getNetworkGridLayout
#'
#' \code{getNetworkGridLayout} return a list x,y coordinates
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param vId the vertex id
#' @author M.Beauvais
#' @export
getNetworkGridLayout<-function(AFMImage, vId) {
  vertexId<-as.numeric(vId)
  y<-floor(vertexId/HASHSIZE)
  x<-vertexId-y*HASHSIZE
  return(data.table(x=x,y=y))
}

#' Does an edge exist ?
#'
#' \code{existsEdge} return TRUE if an edge exists for this vertex id
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param vertexId the vertex id
#' @author M.Beauvais
#' @export
existsEdge<-function(AFMImage, vertexId) {
  # print(vertexId)
  if ((vertexId<1)||(vertexId>(AFMImage@samplesperline+HASHSIZE*(AFMImage@lines-1)))) {
    # print("return FALSE")
    return(FALSE)
  }
  # print(vertexId)
  
  
  coordinates<-getCoordinatesFromVertexId(AFMImage, vertexId)
  # print(coordinnates)
  id<-coordinates[1]+AFMImage@samplesperline*coordinates[2]
  # print(id)
  if (AFMImage@data$h[id]>0) {
    # print("return TRUE")
    return(TRUE)
  }
  # print("return FALSE")
  return(FALSE)
}

#' Get surrounding vertices from x,y coordinates
#'
#' \code{getSurroundingVertexesList} return the vertexId
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param x coordinates in x axis
#' @param y coordinates in y axis
#' @author M.Beauvais
#' @export
getSurroundingVertexesList<-function(AFMImage,x,y) {
  #   print(x)
  #   print(y)
  horizontalWeight<-AFMImage@hscansize/AFMImage@samplesperline
  verticalWeight<-AFMImage@vscansize/AFMImage@lines
  diagWeight<-sqrt((AFMImage@vscansize/AFMImage@lines)^2+(AFMImage@hscansize/AFMImage@samplesperline)^2)
  
  currentVertexId<-getVertexId(AFMImage,x,y)
  vList=data.table()
  #x+1 y
  nearVertexId<-getVertexId(AFMImage,x+1,y) 
  # print(nearVertexId)
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(horizontalWeight)))
  #x+1 y+1
  nearVertexId<-getVertexId(AFMImage,x+1,y+1) 
  #print(existsEdge(AFMImage, nearVertexId))
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(diagWeight)))
  
  #x y+1
  nearVertexId<-getVertexId(AFMImage,x,y+1) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(verticalWeight)))
  
  #x-1 y+1
  nearVertexId<-getVertexId(AFMImage,x-1,y+1) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(diagWeight)))
  
  #x-1 y
  nearVertexId<-getVertexId(AFMImage,x-1,y) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(horizontalWeight)))
  
  #x-1 y-1
  nearVertexId<-getVertexId(AFMImage,x-1,y-1) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(diagWeight)))
  
  #x y-1
  nearVertexId<-getVertexId(AFMImage,x,y-1) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(verticalWeight)))
  
  #x+1 y-1
  nearVertexId<-getVertexId(AFMImage,x+1,y-1) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(diagWeight)))
  return(vList)
}

#' isAdjacentToBetterVertex
#'
#' \code{isAdjacentToBetterVertex} return TRUE if vertex is adjacent to a better vertex
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param x coordinates in x axis
#' @param y coordinates in y axis
#' @author M.Beauvais
#' @export
isAdjacentToBetterVertex<-function(AFMImage,x,y) {
  #   print(x)
  #   print(y)
  
  currentVertexId<-getVertexId(AFMImage,x,y) 
  currentH<-AFMImage@data$h[currentVertexId]
  
  if(currentH<=0) return(FALSE)
  
  #x+1 y
  nearVertexId<-getVertexId(AFMImage,x+1,y) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x+1 y+1
  nearVertexId<-getVertexId(AFMImage,x+1,y+1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x y+1
  nearVertexId<-getVertexId(AFMImage,x,y+1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x-1 y+1
  nearVertexId<-getVertexId(AFMImage,x-1,y+1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x-1 y
  nearVertexId<-getVertexId(AFMImage,x-1,y) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x-1 y-1
  nearVertexId<-getVertexId(AFMImage,x-1,y-1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x y-1
  nearVertexId<-getVertexId(AFMImage,x,y-1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x+1 y-1
  nearVertexId<-getVertexId(AFMImage,x+1,y-1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  return(FALSE)
}

#' gridIgraphPlot
#'
#' \code{gridIgraphPlot} return TRUE if vertex is adjacent to a better vertex
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param g the networks
#' @author M.Beauvais
#' @export
gridIgraphPlot<-function(AFMImage, g){
  # define the layout matrix
  coordinatesVector<-getNetworkGridLayout(AFMImage, V(g)$name)
  #coordinatesVector
  
  l<-matrix(coordinatesVector$x ,byrow = TRUE)
  l<-cbind(l, coordinatesVector$y)
  #l
  
  # plot.igraph(all, layout=All_layout, vertex.size=2, vertex.label=V(All)$name,
  #      vertex.color="green", vertex.frame.color="red", edge.color="grey",  
  #      edge.arrow.size=0.01, rescale=TRUE,vertex.label=NA, vertex.label.dist=0.0,
  #      vertex.label.cex=0.5, add=FALSE,   vertex.label.font=.001)
  plot.igraph(g, layout=l, 
       vertex.shape="circle", vertex.size=2, vertex.label=NA, vertex.color="red", vertex.frame.color="red",
       edge.color="grey"
  )
  
}

#' Calculate iGraph from AFMImage
#'
#' \code{calculateIgraph} return 
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param AFMImageNetworksAnalysis an \code{\link{AFMImageNetworksAnalysis}} from Atomic Force Microscopy
#' @author M.Beauvais
#' @export
calculateIgraph<-function(AFMImage, AFMImageNetworksAnalysis) {
  if (missing(AFMImageNetworksAnalysis)) {
    AFMImageNetworksAnalysis<-NULL
  }
  graphicalUpdate<-FALSE
  graphicalCounter<-0
  
  if (!is.null(AFMImageNetworksAnalysis)&&
      !is.null(AFMImageNetworksAnalysis@updateProgress)&&
      is.function(AFMImageNetworksAnalysis@updateProgress)&&
      !is.null(AFMImageNetworksAnalysis@updateProgress())) {
    graphicalUpdate<-TRUE
    totalLength<-AFMImage@samplesperline*(AFMImage@lines-1)
  }
  
  if (graphicalUpdate) {
    AFMImageNetworksAnalysis@updateProgress(message="1/2 - Generating edges list", value=0)
  }
  print(paste("Generating edge list"))
  
  counter<-1
  #edgeList=data.table()  
  edgeList <- vector("list", AFMImage@samplesperline*AFMImage@lines+1)
  
  for (x in seq(1: AFMImage@samplesperline)) {
    for (y in seq(1: (AFMImage@lines-1))) {
      currentVertexId<-getVertexId(AFMImage,x,y)
      if (existsEdge(AFMImage, currentVertexId)) {
        #edgeList<-rbind(edgeList, getSurroundingVertexesList(AFMImage,x,y))
        edgeList[[counter]] <- getSurroundingVertexesList(AFMImage,x,y)
        counter<-counter+1
      }
      if (graphicalUpdate) {
        graphicalCounter<-graphicalCounter+1
        if (graphicalCounter/100==floor(graphicalCounter/100)) {
        value<-graphicalCounter / totalLength
        text <- paste0(round(graphicalCounter, 2),"/",totalLength)
        AFMImageNetworksAnalysis@updateProgress(value= 0, detail = text)
      }
    }
  }
  }
  
  if (graphicalUpdate) {
    AFMImageNetworksAnalysis@updateProgress(message="2/2 - Generating network", value=0)
  }
  
  newEdgeList<-rbindlist(edgeList)
  el=as.matrix(newEdgeList)
  print(paste("Creating graph"))
  g<-graph_from_edgelist(el[,1:2], directed=FALSE)
  print(paste("Created",counter,"vertices"))
  AFMImageNetworksAnalysis@originalGraph<-g
  return(g)
}

#' getListOfDiameters
#'
#' \code{getListOfDiameters} return 
#' 
#' @param g list of igraph networks
#' @author M.Beauvais
#' @export
getListOfDiameters<-function(g) {
  LIST_OF_DIAMETERS = c()
  listOfGraph=decompose(g)
  for(g in listOfGraph){
    LIST_OF_DIAMETERS=c(LIST_OF_DIAMETERS, diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL))
  }
  return(LIST_OF_DIAMETERS)  
}

#' canBeRemoved
#'
#' \code{canBeRemoved} return 
#' 
#' @param vertexId a vertex id
#' @param g a igraph
#' @param allVertices list of all vertices
#' @param DEGREE_LIMIT_FOR_CANDIDATE_VERTICE degree
#' 
#' @author M.Beauvais
#' @export
canBeRemoved<-function(vertexId, g, allVertices, DEGREE_LIMIT_FOR_CANDIDATE_VERTICE) {
  avList<-adjacent_vertices(g, v=c(vertexId), mode = c("all"))
  avListNew<-unique(avList[[vertexId]]$name)
  found<-NULL
  if (nrow(allVertices[, c("found"):=vertexId %in% avListNew & degree<(DEGREE_LIMIT_FOR_CANDIDATE_VERTICE+1)][found==TRUE])>0) {
    return(FALSE)
  }else{
    return(TRUE)
  }
  
}

#' calculateNetworkSkeleton
#'
#' \code{calculateNetworkSkeleton} return 
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param AFMImageNetworksAnalysis an \code{\link{AFMImageNetworksAnalysis}} from Atomic Force Microscopy
#' @author M.Beauvais
#' @export
calculateNetworkSkeleton<-function(AFMImage, AFMImageNetworksAnalysis) {
  if (missing(AFMImageNetworksAnalysis)) {
    AFMImageNetworksAnalysis<-NULL
    return(new("list"))
  }
  
  g<-AFMImageNetworksAnalysis@originalGraph
  
  graphicalUpdate<-FALSE
  graphicalCounter<-0
  
  if (!is.null(AFMImageNetworksAnalysis)&&
      !is.null(AFMImageNetworksAnalysis@updateProgress)&&
      is.function(AFMImageNetworksAnalysis@updateProgress)&&
      !is.null(AFMImageNetworksAnalysis@updateProgress())) {
    graphicalUpdate<-TRUE
    totalLength<-length(V(g))
    
  }
  
  
  DEGREE_LIMIT_FOR_CANDIDATE_VERTICE=4
  NUMBER_OF_NETWORKS = length(decompose(g))
  LIST_OF_DIAMETERS<-getListOfDiameters(g)
  print(LIST_OF_DIAMETERS)  
  
  
  #   distance_table(g, directed = FALSE)
  #   coreness(g)
  
  
  verticesThatCantBeRemovedList=c()
  print(paste("starting with ", length(V(g)), " vertices"))
  
  if (graphicalUpdate) {
    AFMImageNetworksAnalysis@updateProgress(message="1/1 - removing vertices and edges", value=0)
  }
  
  continueExploration<-TRUE
  while(continueExploration) {
    
    edgeList<-V(g)$name
    
    uniqueVerticesList<-unique(edgeList)
    uniqueVerticesList
    # degree de chaque noeud
    edgeDegreeList<-degree(g, v=uniqueVerticesList, mode = c("all"), loops = FALSE, normalized = FALSE)
    edgeDegreeList
    
    # liste ordonn'e9e croissante des noeuds en fonction du degree
    
    allVertices<-data.table(vertexId=uniqueVerticesList, degree=edgeDegreeList)
    # get-list of adjacent vertices with degree > 2 (can't remove if degree < 2)
    allVertices<-allVertices[order(degree)]
    listOfCandidateVertices<-allVertices[degree>DEGREE_LIMIT_FOR_CANDIDATE_VERTICE]
    
    listOfCandidateVertices<-listOfCandidateVertices[!listOfCandidateVertices$vertexId %in% verticesThatCantBeRemovedList]
    
    continueExploration<-FALSE
    if (nrow(listOfCandidateVertices)>0) {
      
      #             res<-sapply(listOfCandidateVertices$vertexId, canBeRemoved, g=g, allVertices=allVertices, simplify=F)
      #             vMatrix<-as.matrix(res, ncol=2)
      #             
      #             verticesToBeRemoved<-data.table(vertexId= rownames(vMatrix), toBeRemoved= vMatrix[,1])[toBeRemoved==TRUE]$vertexId
      #             print(paste("to be removed",verticesToBeRemoved))
      #             
      #             if (length(verticesToBeRemoved)>0) {
      #               g<-delete_vertices(g, c(verticesToBeRemoved))
      #               #continueExploration<-TRUE
      #               continueExploration<-continueExploration+1
      #             }
      #       
      for (vi in seq(1:nrow(listOfCandidateVertices))){
        onevertexId=listOfCandidateVertices$vertexId[vi]
        if (canBeRemoved(onevertexId, g=g, allVertices=allVertices, DEGREE_LIMIT_FOR_CANDIDATE_VERTICE=DEGREE_LIMIT_FOR_CANDIDATE_VERTICE)) {
          vId<-listOfCandidateVertices$vertexId[vi]
          
          # store the list of adjacent vertices of the node before deleting it
          avList<-unique(adjacent_vertices(g, v=c(vId), mode = c("all"))[[vId]]$name)
          
          
          g<-delete_vertices(g, listOfCandidateVertices$vertexId[vi])
          continueExploration<-TRUE
          
          NEW_LIST_OF_DIAMETERS=getListOfDiameters(g)
          #print(NEW_LIST_OF_DIAMETERS)  
          
          # did the vertex removal split the network or diminish the diameter
          if ((length(decompose(g))>NUMBER_OF_NETWORKS)||(!identical(LIST_OF_DIAMETERS,NEW_LIST_OF_DIAMETERS))) {
            print (paste("should not have removed", vId))
            verticesThatCantBeRemovedList=c(verticesThatCantBeRemovedList, listOfCandidateVertices$vertexId[vi])
            
            g<-g+vertices(as.numeric(vId))
            
            listOfEdges=c()
            for(j in seq(1,length(avList))) {
              listOfEdges=c(listOfEdges, vId, avList[j], avList[j],vId)
            }
            g<-g+edges(listOfEdges)
          }else{
            print("61")
            NEW_LIST_OF_DIAMETERS=getListOfDiameters(g)
            if ((!identical(LIST_OF_DIAMETERS,NEW_LIST_OF_DIAMETERS))) {
              print (paste("should not have removed", vId))
              verticesThatCantBeRemovedList=c(verticesThatCantBeRemovedList, listOfCandidateVertices$vertexId[vi])
              
              g<-g+vertices(as.numeric(vId))
              
              listOfEdges=c()
              for(j in seq(1,length(avList))) {
                listOfEdges=c(listOfEdges, vId, avList[j], avList[j],vId)
          }
          
              g<-g+edges(listOfEdges)
            }
            break
        }
      }
        
      }
      if (graphicalUpdate) {
        graphicalCounter<-graphicalCounter+1
        value<-graphicalCounter / totalLength
        text <- paste0(round(graphicalCounter, 2),"/",totalLength)
        AFMImageNetworksAnalysis@updateProgress(value= 0, detail = text)
      }
      
    }else{
      continueExploration<-FALSE
    }
  }
  print(paste("ending with ", length(V(g)), " vertices"))
  
  AFMImageNetworksAnalysis@skeletonGraph<-g
  
  return(AFMImageNetworksAnalysis)
}

#' Calculate topology image (TBC)
#'
#' \code{getTopologyAFMImage} return the global topological distance
#' 
#' @param BinaryAFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy in a binary format 0 or 1 values for heigths
#' @param AFMImageNetworksAnalysis an \code{\link{AFMImageNetworksAnalysis}} from Atomic Force Microscopy
#' @author M.Beauvais
#' @export
getTopologyAFMImage<-function(BinaryAFMImage, AFMImageNetworksAnalysis){

  filterVector<-unlist(BinaryAFMImage@data$h)
  
  topology<-c()
  
  
  for (x in 1:BinaryAFMImage@samplesperline) {
    for (y in 1:BinaryAFMImage@lines) {
      if(x==1) {
        bX=seq(from=0, to=BinaryAFMImage@samplesperline-1, by=1)
      }else{
        if (x==BinaryAFMImage@samplesperline) {
          bX=seq(from=x-1, to=0, by=-1)
        }else{
          bX=seq(from=x-1, to=0, by=-1)
          bX=c(bX, seq(from=1, to=BinaryAFMImage@samplesperline-x, by=1))
        }
      }
      # bX
      
      if(y==1) {
        bY=seq(from=0, to=BinaryAFMImage@lines-1, by=1)
      }else{
        if (y==BinaryAFMImage@lines) {
          bY=seq(from=y-1, to=0, by=-1)
        }else{
          bY=seq(from=y-1, to=0, by=-1)
          bY=c(bY, seq(from=1, to=BinaryAFMImage@lines-y, by=1))
        }
      }
      # bY
      
      
      bX=BinaryAFMImage@hscansize*bX
      bY=BinaryAFMImage@vscansize*bY
      
      bX<-matrix(rep(bX,BinaryAFMImage@lines), ncol=BinaryAFMImage@lines, byrow=TRUE )
      bY<-matrix(rep(bY,BinaryAFMImage@samplesperline), ncol=BinaryAFMImage@samplesperline, byrow=FALSE )
      
      nm=as.numeric(1/sqrt(bX^2+bY^2))
      nm[is.infinite(nm)]<-0
      #nm*filterVector
      res<-sum(nm*filterVector)
      topology<-c(topology,res)
      #print(res)
      
    }
  }
  
  
  scanby<-BinaryAFMImage@scansize/BinaryAFMImage@samplesperline
  endScan<-BinaryAFMImage@scansize*(1-1/BinaryAFMImage@samplesperline)
  
  topologyAFMImage<-AFMImage(
    data = data.table(x = rep(seq(0,endScan, by= scanby), times = BinaryAFMImage@lines),
                      y = rep(seq(0,endScan, by= scanby), each = BinaryAFMImage@samplesperline),
                      h = topology),
    samplesperline = BinaryAFMImage@samplesperline, lines = BinaryAFMImage@lines,
    vscansize = BinaryAFMImage@vscansize, hscansize = BinaryAFMImage@hscansize, scansize = BinaryAFMImage@scansize,
    fullfilename = BinaryAFMImage@fullfilename )
  
  
  
  return(topologyAFMImage)
  
}

#' get a segment of points thanks to Bresenham line algorithm
#'
#' \code{getBresenham2DSegment} return the Bresenham segment in 2D from extremities coordinates
#' 
#' @param x1 abscissa coordinates of the first point
#' @param y1 ordinate coordinates of the first point
#' @param x2 abscissa coordinates of the second point
#' @param y2 ordinate coordinates of the second point
#' 
#' @author M.Beauvais
#' @export

getBresenham2DSegment<-function(x1, y1, x2, y2) {
  resX=c()
  resY=c()
  
  dx<-x2-x1
  dy<-y2-y1
  
  #print(paste("getBresenham2DSegment",dx,dy))
  
  if (dx !=0) {
    if (dx > 0) {
      if (dy !=0) {
        if (dy > 0) {
          if (dx >= dy) {
            e<-dx
            dx <- e  * 2 
            dy <- dy * 2  
            while(TRUE){  
              resX=c(resX,x1); resY=c(resY, y1)
              x1 <- x1 + 1
              if (x1 == x2) break
              e <- e - dy
              if (e < 0) {
                y1 <- y1 + 1
                e <- e + dx 
              }
            }
          } else {
            e <- dy
            dy <- e * 2
            dx <- dx * 2 
            while(TRUE){ 
              resX=c(resX,x1); resY=c(resY, y1)
              y1 <- y1 + 1
              if (y1 == y2) break
              e <- e - dx
              if (e < 0) {
                x1 <- x1 + 1 
                e <- e + dy
              }
            }
          }
        }else if (dy < 0){ # dy < 0 (et dx > 0)
          
          
          if (dx >= -dy) {
            
            e <- dx
            dx <- e * 2
            dy <- dy * 2
            
            while(TRUE){  
              resX=c(resX,x1); resY=c(resY, y1)
              x1 <- x1 + 1
              if (x1 == x2) break
              e <- e + dy
              if (e < 0) {
                y1 <- y1 - 1 
                e <- e + dx
              }
            }
          } else{
            
            e <- dy
            dy <- e * 2 
            dx <- dx * 2
            while(TRUE){  
              resX=c(resX,x1); resY=c(resY, y1)
              y1 <- y1 - 1
              if (y1 == y2) break
              e <- e - dx
              if (e > 0) {
                x1 <- x1 + 1
                e <- e - dy
              }
            }
          }
          
        }
      }  else if (dy == 0){ # dy = 0 (et dx > 0)
        while(x1 != x2) {
          resX=c(resX,x1); resY=c(resY, y1) 
          x1 <- x1 + 1
        }
      }
    }else if (dx<0) {  # dx < 0
      dy <- y2 - y1
      if (dy != 0) {
        if (dy > 0) {
          if (-dx >= dy) {
            e <- dx
            dx <- e * 2 
            dy <- dy * 2  
            while(TRUE){
              resX=c(resX,x1); resY=c(resY, y1) 
              x1 <- x1 - 1
              if (x1 == x2) break
              e <- e + dy
              if (e >= 0) {
                y1 <- y1 + 1 
                e <- e + dx 
              }
            }
          }else{
            e <- dy
            dy <- e * 2
            dx <- dx * 2 
            while(TRUE){ 
              resX=c(resX,x1); resY=c(resY, y1) 
              y1 <- y1 + 1
              if ( y1 == y2) break 
              e <- e + dx
              if (e <= 0) {
                x1 <- x1 - 1  
                e <- e + dy 
              }
            }
          }
        }else if(dy <0) {  # dy < 0 (et dx < 0)
          if (dx <= dy) {
            e <- dx
            dx <- e * 2 
            dy <- dy * 2  
            while(TRUE){  
              resX=c(resX,x1); resY=c(resY, y1)
              x1 <- x1 - 1
              if (x1 == x2) break
              e <- e - dy
              if (e >= 0) {
                y1 <- y1 - 1
                e <- e + dx 
              }
            }
          } else { 
            e <- dy
            dy <- e * 2 
            dx <- dx * 2 
            
            while(TRUE){
              resX=c(resX,x1); resY=c(resY, y1)
              y1 <- y1 - 1
              if ( y1 == y2 ) break
              e <- e - dx
              if (e >= 0) {
                x1 <- x1 - 1
                e <- e + dy
              }
            }
          }
        } 
      } else if (dy==0) {  # dy = 0 (et dx < 0)
        while(x1!=x2) {
          resX=c(resX,x1); resY=c(resY, y1)
          x1 <- x1 - 1
        }
      }
    }
  } else if (dx==0) {  # dx = 0
    dy <- y2 - y1
    if (dy != 0) {
      if (dy > 0) {
        while(y1 != y2) {
          resX=c(resX,x1); resY=c(resY, y1)
          y1 <- y1 + 1
        } 
        
      } else if (dy < 0) { # dy < 0 (et dx = 0)
        while(y1!=y2) {
          resX=c(resX,x1); resY=c(resY, y1)
          y1 <- y1 - 1
        }
        
      }
      
    }
    
  }
  resX=c(resX,x2); resY=c(resY, y2)
  pts = data.table(x=resX, y=resY)
  
  return(pts)
}
