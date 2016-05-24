triangulate_native <- function(P, PB, PA, S, SB,H, TR, flags) {
  ## It is necessary to check for NAs and NaNs, as the triangulate C
  ## code crashes if fed with them
  
  P  <- as.matrix(P)
  PB <- as.integer(PB)
  PA <- as.matrix(PA)
  S  <- as.matrix(S)
  SB <- as.integer(SB)
  H  <- as.matrix(H)
  TR  <- as.matrix(TR)
  
  storage.mode(P)  <- "double"
  storage.mode(PA) <- "double"
  #storage.mode(PB) <- "integer"
  storage.mode(S)  <- "integer"
  #storage.mode(SB) <- "integer"
  storage.mode(H)  <- "double"
  storage.mode(TR) <- "integer"
  storage.mode(flags) <- 'character'
  ## Call the main routine
  out <- .Call("R_triangulate_native",
               t(P),
               PB,
               PA,
               t(S),
               SB,
               H,
               t(TR),
               flags,
               PACKAGE="fdaPDE")
  names(out) <- c("P", "PB", "PA", "T", "S", "SB", "E", "EB","TN", "VP", "VE", "VN", "VA")
  class(out) <- "triangulation"
  return(out)
}

#' Create a triangular mesh
#' 
#' @param nodes A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.
#' @param nodesattributes A matrix with #nodes rows containing nodes' attributes. 
#' These are passed unchanged to the output. If a node is added during the triangulation process or mesh refinement, its attributes are computed  
#' by linear interpolation using the attributes of neighboring nodes. This functionality is for instance used to compute the value 
#' of a Dirichlet boundary condition at boundary nodes added during the triangulation process.
#' @param segments A #segments-by-2 matrix. Each row contains the row's indices in \code{nodes} of the vertices where the segment starts from and ends to.
#' Segments are edges that are not splitted during the triangulation process. These are for instance used to define the boundaries
#' of the domain. If this is input is NULL, it generates a triangulation over the
#' convex hull of the points specified in \code{nodes}.
#' @param holes A #holes-by-2 matrix containing the x and y coordinates of a point internal to each hole of the mesh. These points are used to carve holes
#' in the triangulation, when the domain has holes.
#' @param triangles A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix.
#' This option is used when a triangulation is already available. It specifies the triangles giving the row's indices in \code{nodes} of the triangles' vertices and (when \code{nodes} = 2) also if the triangles' edges midpoints. The triangles' vertices and midpoints are ordered as described 
#' at \cr https://www.cs.cmu.edu/~quake/triangle.highorder.html.
#' In this case the function \code{create.MESH.2D} is used to produce a complete MESH2D object. 
#' @param order Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints). 
#' These are
#' respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements. Default is \code{order} = 1.
#' @param verbosity This can be '0', '1' or '2'. It indicates the level of verbosity in the triangulation process. When \code{verbosity} = 0 no message is returned
#' during the triangulation. When \code{verbosity} = 2 the triangulation process is described step by step by displayed messages.
#' Default is \code{verbosity} = 0.
#' @description This function is a wrapper of the Triangle library (http://www.cs.cmu.edu/~quake/triangle.html). It can be used
#' to create a triangulation of the domain of interest starting from a list of points, to be used as triangles' vertices, and a list of segments, that define the domain boundary. The resulting
#' mesh is a Constrained Delaunay triangulation. This is constructed in a way to preserve segments provided in the input \code{segments} without splitting them. This imput can be used to define the boundaries
#' of the domain. If this imput is NULL, it generates a triangulation over the
#' convex hull of the points.
#' @usage create.MESH.2D(nodes, nodesattributes = NA, segments = NA, holes = NA, 
#'                      triangles = NA, order = 1, verbosity = 0)
#' @seealso \code{\link{refine.MESH.2D}}, \code{\link{create.FEM.basis}}
#' @return An object of the class MESH2D with the following output:
#' \item{\code{nodes}}{A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.}
#' \item{\code{nodesmarkers}}{A vector of length #nodes, with entries either '1' or '0'. An entry '1' indicates that the corresponding node is a boundary node; an entry '0' indicates that the corresponding node is not a boundary node.}
#' \item{\code{nodesattributes}}{nodesattributes A matrix with #nodes rows containing nodes' attributes. 
#' These are passed unchanged to the output. If a node is added during the triangulation process or mesh refinement, its attributes are computed  
#' by linear interpolation using the attributes of neighboring nodes. This functionality is for instance used to compute the value 
#' of a Dirichlet boundary condition at boundary nodes added during the triangulation process.}
#' \item{\code{triangles}}{A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix.
#' This option is used when a triangulation is already available. It specifies the triangles giving the indices in \code{nodes} of the triangles' vertices and (when \code{nodes} = 2) also if the triangles' edges midpoints. The triangles' vertices and midpoints are ordered as described 
#' at  \cr https://www.cs.cmu.edu/~quake/triangle.highorder.html.}
#' \item{\code{segmentsmarker}}{A vector of length #segments with entries either '1' or '0'. An entry '1' indicates that the corresponding element in \code{segments} is a boundary segment;  
#' an entry '0' indicates that the corresponding segment is not a boundary segment.}
#' \item{\code{edges}}{A #edges-by-2 matrix containing all the edges of the triangles in the output triangulation. Each row contains the row's indices in \code{nodes}, indicating the nodes where the edge starts from and ends to.}
#' \item{\code{edgesmarkers}}{A vector of lenght #edges with entries either '1' or '0'. An entry '1' indicates that the corresponding element in \code{edge} is a boundary edge;  
#' an entry '0' indicates that the corresponding edge is not a boundary edge.}
#' \item{\code{neighbors}}{A #triangles-by-3 matrix. Each row contains the indices of the three neighbouring triangles. An entry '-1' indicates that 
#' one edge of the triangle is a boundary edge.}
#' \item{\code{holes}}{A #holes-by-2 matrix containing the x and y coordinates of a point internal to each hole of the mesh. These points are used to carve holes
#' in the triangulation, when the domain has holes.}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints). 
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements. Default is \code{order} = 1.}
#' @examples 
#' ## Upload the Meuse data
#' data(MeuseData)
#' ## Create a triangulation on the convex hull of these data,
#' ## where each data location is a triangle vertex
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], order = 1)
#' ## Plot the mesh
#' plot(mesh)
#' ## Upload a domain boundary for these data
#' data(MeuseBorder)
#' ## Create a constrained Delaunay triangulation with the provided boundary 
#' ## where each datalocation is a triangle vertex
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = 1)
#' ## Plot the mesh
#' plot(mesh)

create.MESH.2D <- function(nodes, nodesattributes = NA, segments = NA, holes = NA, triangles = NA, order = 1, verbosity = 0)
{ 
  ##########################
  ###   Input checking   ###
  ##########################
  
  # Triangle finds out which are on the border (see https://www.cs.cmu.edu/~quake/triangle.help.html)
  nodesmarkers = vector(mode = "integer", 0)
  segmentsmarkers = vector(mode = "integer", 0)
  
  nodes = as.matrix(nodes)
  if (ncol(nodes) != 2)
    stop("Matrix of nodes should have 2 columns")
  if (anyDuplicated(nodes))
    stop("Duplicated nodes")
  
  ## If attributes not specified, set them to a matrix with zero columns
  if (any(is.na(nodesattributes))) {
    nodesattributes <- matrix(0, nrow(nodes), 0)
  }else{
    nodesattributes <- as.matrix(nodesattributes)
    if (nrow(nodesattributes) != nrow(nodes))
      stop("Point attribute matrix \'nodesattributes\' does not have same number of rows the point matrix \'nodes\'")
  }
  
  ## If boundary nodes not specified, set them to 0
#   if (any(is.na(nodesmarkers))) {
#     nodesmarkers <- vector(mode = "integer", 0)
#   }else{
#     nodesmarkers = as.vector(nodesmarkers)
#   }
    
  ## Deal with segments
  if (any(is.na(segments))) {
    segments <- matrix(0, 0, 2)
  } else {
    segments <- as.matrix(segments)
    if (ncol(segments) != 2) {
      stop("Matrix of segments should have 2 columns")
    }
  }
  
  ## If boundary segments not specified, set them to 0
#   if (any(is.na(segmentsmarkers))) {
#     segmentsmarkers <- vector(mode = "integer", 0)
#   }else{
#     segmentsmarkers = as.vector(segmentsmarkers)
#   }
  
  ## If hole not specified, set it to empty matrix
  if (any(is.na(holes)))
    holes <- matrix(0, 0, 2)
  holes = as.matrix(holes)
  
  ## If triangles are not already specified
  if(any(is.na(triangles)))
    triangles = matrix(0,nrow = 0, ncol = 3)
  triangles = as.matrix(triangles)
  
  ## Set meshing parameters ##
  flags="ven"
  if(nrow(segments) == 0){
    flags = paste(flags,"c",sep = '')
  }
  
  if(nrow(segments)>0){
    flags = paste(flags,"p",sep = '')
  }
  
  #If order=2 add flag for second order nodes
  if(order == 2){
    flags = paste(flags,"o2",sep = '')
  }
  if(order < 1 || order >2){
    print('Order must be 1 or 2')
  }
  
  if(nrow(triangles) > 0){
    flags = paste(flags,"r",sep = '')
  }
  
  if (verbosity == 0) {
    flags = paste(flags,"Q",sep = '')
  }
  if (verbosity == 1) {
    flags = paste(flags,"V",sep = '')
  }
  if (verbosity == 2) {
    flags = paste(flags,"VV",sep = '')
  }
  
  out<-NULL
  #If triangles is null it makes the trianglulation
  #If triangle is not null it makes a refinement with no parameter, to compose the mesh object
  out <- triangulate_native(  
    nodes,
    nodesmarkers,
    nodesattributes,
    segments,
    segmentsmarkers,
    t(holes),
    triangles,
    flags
  )
  
  names(out)[1]<-"nodes"
  names(out)[2]<-"nodesmarkers"
  names(out)[3]<-"nodesattributes"
  names(out)[4]<-"triangles"
  names(out)[5]<-"segments"
  names(out)[6]<-"segmentsmarkers"
  names(out)[7]<-"edges"
  names(out)[8]<-"edgesmarkers"
  names(out)[9]<-"neighbors"
  
  out[13]<-NULL
  out[12]<-NULL
  out[11]<-NULL
  out[10]<-NULL
  
  out[[10]] = holes
  names(out)[10]<-"holes"
  out[[11]] = order
  names(out)[11]<-"order"
  
  class(out)<-"MESH2D"
  
  return(out)
}

#' Refine a triangular mesh
#' 
#' @param mesh A MESH2D object representing the triangular mesh, created by \link{create.MESH.2D}.
#' @param minimum_angle A scalar specifying a minimun value for the triangles angles.
#' @param maximum_area A scalar specifying a maximum value for the triangles areas.
#' @param delaunay A boolean parameter indicating whether or not the output mesh should satisfy the Delaunay condition.
#' @param verbosity This can be '0', '1' or '2'. It indicates the level of verbosity in the triangulation process.
#' @description This function refines a Constrained Delaunay triangulation into a Conforming Delaunay triangulation. This is a wrapper of the Triangle library (http://www.cs.cmu.edu/~quake/triangle.html). It can be used to 
#' refine a mesh created previously with \link{create.MESH.2D}. The algorithm can add Steiner points (points through which the \code{segments} are splitted)
#' in order to meet the imposed refinement conditions.
#' @usage refine.MESH.2D(mesh, minimum_angle, maximum_area, delaunay, verbosity)
#' @seealso \code{\link{create.MESH.2D}}, \code{\link{create.FEM.basis}}
#' @return A MESH2D object representing the refined triangular mesh,  with the following output:
#' \item{\code{nodes}}{A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.}
#' \item{\code{nodesmarkers}}{A vector of length #nodes, with entries either '1' or '0'. An entry '1' indicates that the corresponding node is a boundary node; an entry '0' indicates that the corresponding node is not a boundary node.}
#' \item{\code{nodesattributes}}{nodesattributes A matrix with #nodes rows containing nodes' attributes. 
#' These are passed unchanged to the output. If a node is added during the triangulation process or mesh refinement, its attributes are computed  
#' by linear interpolation using the attributes of neighboring nodes. This functionality is for instance used to compute the value 
#' of a Dirichlet boundary condition at boundary nodes added during the triangulation process.}
#' \item{\code{triangles}}{A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix.
#' This option is used when a triangulation is already available. It specifies the triangles giving the row's indices in \code{nodes} of the triangles' vertices and (when \code{nodes} = 2) also if the triangles' edges midpoints. The triangles' vertices and midpoints are ordered as described 
#' at \cr  https://www.cs.cmu.edu/~quake/triangle.highorder.html.}
#' \item{\code{edges}}{A #edges-by-2 matrix. Each row contains the row's indices of the nodes where the edge starts from and ends to.}
#' \item{\code{edgesmarkers}}{A vector of lenght #edges with entries either '1' or '0'. An entry '1' indicates that the corresponding element in \code{edge} is a boundary edge;  
#' an entry '0' indicates that the corresponding edge is not a boundary edge.}
#' \item{\code{neighbors}}{A #triangles-by-3 matrix. Each row contains the indices of the three neighbouring triangles. An entry '-1' indicates that 
#' one edge of the triangle is a boundary edge.}
#' \item{\code{holes}}{A #holes-by-2 matrix containing the x and y coordinates of a point internal to each hole of the mesh. These points are used to carve holes
#' in the triangulation, when the domain has holes.}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints). 
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements. Default is \code{order} = 1.}
#' @examples 
#' ## Upload the Meuse data and a domain boundary for these data
#' data(MeuseData)
#' data(MeuseBorder)
#' ## Create a Constrained Delaunay triangulation
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = 1)
#' ## Plot the mesh
#' plot(mesh)
#' ## Refine the triangulation
#' mesh_refine <- refine.MESH.2D(mesh, minimum_angle = 30, maximum_area = 10000)
#' plot(mesh_refine)

refine.MESH.2D<-function(mesh, minimum_angle = NA, maximum_area = NA, delaunay = FALSE, verbosity = 0)
{ 
  flags="rpven" 
  
  if(!is.na(minimum_angle)){
    flags <- paste(flags, "q", sprintf("%.12f", minimum_angle), sep='')
  }
  
  if(!is.na(maximum_area)){
    flags <- paste(flags, "a", sprintf("%.12f", maximum_area), sep='')
  }
  
  if(delaunay){
    flags <- paste(flags, "D", sep='')
  }
  
  if(mesh$order==2){
    flags <- paste(flags, "o2", sep='')
  }
  
  if (verbosity == 0) {
    flags = paste(flags,"Q",sep = '')
  }
  if (verbosity == 1) {
    flags = paste(flags,"V",sep = '')
  }
  if (verbosity == 2) {
    flags = paste(flags,"VV",sep = '')
  }
  
  # Triangle finds out which are on the border (see https://www.cs.cmu.edu/~quake/triangle.help.html)
  mesh$nodesmarkers = vector(mode = "integer", 0)
  mesh$segmentsmarkers = vector(mode = "integer", 0)
  
  out<-NULL
  #If triangles is null it makes the trianglulation
  #If triangle is not null it makes a refinement with no parameter, to compose the mesh object
  out <- triangulate_native(
    mesh$nodes,
    mesh$nodesmarkers,
    mesh$nodesattributes,
    mesh$segments,
    mesh$segmentmarkers,
    t(mesh$holes),
    mesh$triangles,
    flags
  )
  
  names(out)[1]<-"nodes"
  names(out)[2]<-"nodesmarkers"
  names(out)[3]<-"nodesattributes"
  names(out)[4]<-"triangles"
  names(out)[5]<-"segments"
  names(out)[6]<-"segmentsmarkers"
  names(out)[7]<-"edges"
  names(out)[8]<-"edgesmarkers"
  names(out)[9]<-"neighbors"
  
  out[13]<-NULL
  out[12]<-NULL
  out[11]<-NULL
  out[10]<-NULL
  
  out[[10]] = mesh$holes
  names(out)[10]<-"holes"
  out[[11]] = mesh$order
  names(out)[11]<-"order"
  
  class(out)<-"MESH2D"
  
  return(out)
}

#' Plot a MESH2D object
#' 
#' @param x A MESH2D object defining the triangular mesh, as generated by \code{create.Mesh.2D} or \code{refine.Mesh.2D}.
#' @param ... Arguments representing graphical options to be passed to \link[graphics]{par}.
#' @description Plot a mesh MESH2D object, generated by \code{create.MESH.2D} or \code{refine.MESH.2D}. Circles indicate the mesh nodes.
#' @usage \method{plot}{MESH2D}(x, ...)
#' @examples 
#' ## Upload the Meuse data and a domain boundary
#' data(MeuseData)
#' data(MeuseBorder)
#' ## Create a triangular mesh with the provided boundary
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = 1)
#' ## Plot it
#' plot(mesh)
plot.MESH2D<-function(x, ...)
{
  plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
  segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
           x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
  segments(x$nodes[x$segments[,1],1], x$nodes[x$segments[,1],2],
           x$nodes[x$segments[,2],1], x$nodes[x$segments[,2],2], col="red", ...)
}
