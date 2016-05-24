#  printing and graphing functions for objects of class mvmesh
#######################################################################
print.mvmesh <- function( x, ... ) {
# quick method to print out contents of mesh x

if( class(x) != "mvmesh") warning( "x is not a mesh object" )

print(str(x)) }
#######################################################################
plot.mvmesh <- function( x, new.plot=TRUE, show.points=FALSE, show.edges=TRUE, show.faces=FALSE, 
                         show.labels=FALSE, label.values=NULL, ...) {
#  Plot the mesh x in n=2  or n=3 dimensions
#  x is a list of class "mvmesh"
#  If new.plot=TRUE, a new plot is started; if FALSE, graphing is added to an existing plot
#  If show.faces=TRUE, faces of the simplex are shown; otherwise edges are shown
#  If show.points=TRUE, vertices of the simplex are shown
#  If show.labels=TRUE, the simplices are labeled on the plot
#  label.values are the labels to show; defaults to 1,2,3,... if not specified and show.label=TRUE
#  ... optional arguments to plot, e.g. col, type, lwd, etc.

# check input & define values
if( class(x) != "mvmesh" ) stop( "x is not a mesh object" )
n <- ncol(x$S)
if( (n < 2) | (n > 3) ) stop( "x is not a 2 or 3 dimensional object" )
stopifnot( nrow(x$S) >= 1 )
nS <- dim(x$S)[3]   
if( is.null(label.values) ) { label.values <- 1:nS }
stopifnot( length(label.values)==nS )

# plot the mesh
if (n==2) {
  if (new.plot) {
    # initialize the plot region to show all points
    r1 <- range(as.double(x$S[,1,]))
    r2 <- range(as.double(x$S[,2,]))
    plot( r1, r2, xlab="", ylab="", type="n",...)  
  }
  if (show.points) {
    points( x$V[,1], x$V[,2] , ... )
  }
  for (k in 1:nS ) {
    DrawSimplex2d( x$S[,,k], label.values[k], show.labels, x$mvmesh.type, show.edges, show.faces, ...)
  }
} else { # n==3
   if (new.plot) {
    open3d( )
  }  
  old.par <- par3d( skipRedraw=TRUE ) # speed up plot
  on.exit( par3d( old.par ) ) # restore old settings when done 
  if (show.points) {
    points3d( x$V[,1], x$V[,2], x$V[,3] , ... )
  }
  for (k in 1:nS ) {
    DrawSimplex3d( x$S[,,k], label.values[k], show.labels, x$mvmesh.type, show.edges, show.faces,  ...)
  }
}
}
#####################################################################
DrawSimplex2d <- function( S, label, show.labels, mvmesh.type, show.edges=TRUE, show.faces=FALSE, ...){
#  Draw a single 2d simplex; if show.label=TRUE, add a label 'label' at the middle of the simplex

if ( mvmesh.type==7L ) {
  # for a rectangular mesh, draw a rectangle; need to reorder vertices
  # warning: this reordering depends on the order in which MeshRectangular defines vertices!
  jj <- c(1,2,4,3,1) 
} else {
  jj <- c(1:nrow(S),1)
}      

if (show.edges) {
  lines(S[jj,1], S[jj,2], ... )
}

if (show.faces) {
  polygon( S[jj,1], S[jj,2], ... )
}

if (show.labels) {
  # label simplices
  mid <- colMeans(S)
  text( mid[1], mid[2], label, ... )
}

}
#######################################################################
DrawSimplex3d <- function( S, label, show.labels, mvmesh.type, show.edges=TRUE, show.faces=FALSE, ...){
#  Draw a single 3d simplex.  Default is a wireframe, e.g. lines connecting the vertices; 
#  set show.faces=TRUE to show solid faces.  If show.labels=TRUE, add a 
#  label 'num' at the middle of the simplex

if (mvmesh.type==7L) {
  # rectangular mesh has to be drawn differently
  face.grouping <- matrix( c(1,2,6,5,1,  2,4,8,6,2,  4,3,7,8,4,  
                             3,7,5,1,3,  5,6,8,7,5,  1,2,4,3,1),
                             byrow=TRUE,nrow=6,ncol=5)
  for (k in 1:6) {
    # draw face k
    jj <- face.grouping[k,]
    if (show.edges) { lines3d( S[jj,1],S[jj,2],S[jj,3], ... ) }
    if (show.faces) { jj <- jj[1:4]; rgl.quads( S[jj,1], S[jj,2], S[jj,3], ... ) }
  }
} else {
  if(mvmesh.type==9L) {
    # polar sphere
    jj <- c(1,2,4,3,1)
    if (show.edges) { lines3d( S[jj,1], S[jj,2], S[jj,3], ... ) }
    if (show.faces) { warning("show.faces does not work in this case") }
  } else {
    m <- nrow(S)  
    if (show.edges) { 
      for (i in 1:(m-1)) {
        for (j in (i+1):m) {
          lines3d(S[c(i, j),1],S[c(i, j),2],S[c(i, j),3],...) 
        }  
      }
    }
    if (show.faces) { 
      if (nrow(S)==3) { rgl.triangles( S, ... ) } 
      if (nrow(S)==4) { 
        face.grouping <- matrix( c(1,2,3,  1,2,4,  1,3,4,  2,3,4), byrow=TRUE, nrow=4, ncol=3)
        for (k in 1:4) { 
          jj <- face.grouping[k, ]
          rgl.triangles( S[jj,], ... ) 
        }
      } 
      if (nrow(S) > 4) { warning("show.faces does not work in this case") }         
    }
  }
}

if (show.labels) {
  # label simplices
  mid <- colMeans(S)
  text3d( mid, texts=label, ... )
}
}

