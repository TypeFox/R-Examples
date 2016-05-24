#  R package to generate meshes and grids in n-dimensional Euclidean space.
#
#  Initial functions written in May 2011 by John Nolan (jpnolan@american.edu) 
#  Modified several times thru March 2015 by John Nolan
#     Many new functions added, some bugs fixed.
#
#  This research was supported by an agreement with Cornell University, Operations 
#  Research & Information Engineering, under contract W911NF-12-1-0385 from the Army 
#  Research Development and Engineering Command.
#
#  These R functions provide methods of subdividing and manipulating
#  grids (a list of points) and meshes (a list of points with grouping 
#  information telling which points are the vertices of a simplex). 
#
#####################################################################
UnitSimplex <- function( n, k=1 ) {
#   Generate a uniform mesh on the unit simplex in R^n by a
#   k-subdivision of each edge.  The unit simplex is the
#   (n-1) dimensional surface given by the convex
#   hull of the standard basis vectors e[1],...,e[n].
#
# Output is a list of class "mvmesh"

n <- as.integer( n )
k <- as.integer( k )
# calculate the subdivision for a (n-1)-dim. simplex
T <- EdgeSubdivision( n=n-1, k=k )

# define the original simplex: vertices at (1,0,0,...),(0,1,0,...), etc.
S <- diag( rep(1.0,n) )  

# calculate coordinates of the sub simplices
subS <- array( 0.0, dim=c(n,n,k^(n-1)) )
if (k==1) { subS[,,1] <- S }  # special case of no subdivision
else {
  for (i in 1:k^(n-1)) {  
    subS[,,i] <- SimplexCoord(S,T[,,i])
  }
}

# calculate the (unique) vertices among all sub simplices
a <- SVIFromColor(S,T)

mesh <- list(type="UnitSimplex",mvmesh.type=1L,n=n,m=n-1L,vps=n,S=subS,V=a$V,SVI=a$SVI,k=k)
class(mesh) <- "mvmesh"
return(mesh) }
#####################################################################
SolidSimplex <- function( n, k=1 ) {
#   Generate a uniform mesh on the canonical simplex in R^n by a
#   k-subdivision of each face.  The canonical simplex is an
#   n dimensional solid; it is the convex hull of the standard basis 
#   vectors e[1],...,e[n] and the origin.
#
# Output is a list of class "mvmesh"

n <- as.integer( n )
k <- as.integer( k )
# calculate the subdivision for n-dim. simplex
T <- EdgeSubdivision( n=n, k=k )
M <- dim(T)[3]

# define the original simplex: vertices at (1,0,0,...),(0,1,0,...), etc.
S <- rbind( diag( rep(1,n) ), rep(0,n) )  # size (n+1) x n

# calculate coordinates of the sub-simplices
subS <- array( 0.0, c(n+1,n,M) )
if (k==1) { subS[,,1] <- S }  # special case of no subdivision
else {
  for (i in 1:M) {
    subS[,,i] <- SimplexCoord(S,T[,,i])
  }
}

# calc ulate the (unique) vertices among all sub simplices
a <- SVIFromColor(S,T)

mesh <- list(type="SolidSimplex",mvmesh.type=2L,n=n,m=n,vps=n+1L,S=subS,V=a$V,SVI=a$SVI,k=k)
class(mesh) <- "mvmesh"
return(mesh) }
#####################################################################
UnitSphere <- function( n, k=1, method="dyadic", p=2, positive.only=FALSE ) {
# Generate an "approximately uniform" mesh on the l^p unit ball in R^n. 
# This mesh is exactly uniform when p=1; in other cases, the curvature of 
# the ball makes spacing unequal. 
#
# Input values:
#   n = dimension of the space
#   method = "edgewise" or "dyadic".  
#        "edgewise" uses edgewise subdivision of the unit simplex, and projects onto the sphere. 
#        "dyadic" uses recursive subdivisions by 2 to get a more uniform mesh
#   k = how many times to subdivide.  In the "edgewise' case, this is the number 
#          of subdivisions of unit simplex before projecting.  In the "dyadic" case, 
#          it is the number of recursive subdivisions.
#   p = power used to compute l-p norm: |(x[1],...,x[n])|= (sum(abs(x)^p)^(1/p)
#   positive.only = TRUE to restrict to first quadrant/octant/orthant, = FALSE gives whole ball
#
# Output is a list of class "mvmesh"
#

n <- as.integer( n )
k <- as.integer( k )
stopifnot( n > 1, k >= 0, method %in% c("dyadic","edgewise"), length(p)==1, 
           p > 0, is.logical(positive.only) )

if (method=="edgewise") {
  a <- UnitSphereEdgewise( n=n, k=k, p=p, positive.only=positive.only )
} else {
  a <- UnitSphereDyadic(   n=n, k=k, p=p, positive.only=positive.only ) 
}
return(a) }
#####################################################################
UnitSphereEdgewise <- function( n, k, p, positive.only ) {
# Generate an "approximately uniform" mesh on the l^p unit ball in R^n, 
# based on an order k edge subdivision.  

n <- as.integer( n )
k <- as.integer( k )
stopifnot(k > 0)
#  compute the part of the ball in the positive orthant
a <- UnitSimplex( n=n, k=k )
V <- a$V
M <- nrow(V)
SVI <- a$SVI  

# scale to be on the unit ball in $\ell^p$
if (p != 1.0) { 
  r <- LpNorm( V, p )
  for (i in 1:M) {
    V[i,] <- V[i,]/r[i]
  }
}

# if !positive.only, rotate & reflect first octant to cover all other
# octants.  Check for duplicate vertices and edges along axes
if (!positive.only ) {
  # reflect first orthant into all other orthants
  newptr <- rep( 0L, M )     
  for (i in 1:((2^n)-1)) {
    y <- ConvertBase( i, 2, n )
    signs <- (-1)^y

    # first compute new vertices (some on axis will repeat)
    newM <- nrow(V)
    newV <- rbind( V, matrix(NA,nrow=M, ncol=n ) )  
    for (j in 1:M) {
      tmp <- signs*V[j,]
      l <- MatchRow( tmp, newV )    
      if ( length(l) == 0 ) {
        # tmp is new vector; add to matrix
        newM <- newM + 1
        newV[newM,] <- tmp
        newptr[j] <- newM 
      } else {             
        # tmp is already in newV, just save pointer     
        newptr[j] <- l[1]
      }
    }
    V <- newV[1:newM,] # truncate unused part of V matrix

    # compute simplices in the new orthant
    newSVI <- matrix( 0L, nrow=n, ncol=ncol(a$SVI) )
    for (j in 1:ncol(newSVI)) {
      for (m in 1:n) {
        newSVI[m,j] <- newptr[ SVI[m,j] ]
      }
    }
    SVI <- cbind(SVI,newSVI)    
  }
}

# calculate coordinates of the sub simplices
nSVI <- ncol(SVI)
subS <- array( 0.0, c(n,n,nSVI) )
for (i in 1:nSVI) { 
  for (j in 1:n) {   
    subS[j,,i] <- V[ SVI[j,i], ]
  }
}

mesh <- list( type="UnitSphere, edgewise",mvmesh.type=3L,n=n,m=n-1L,vps=n,  
       S=subS, V=V, SVI=SVI, k=k, p=p , positive.only=positive.only)
class(mesh) <- "mvmesh"
return( mesh ) }
#####################################################################
UnitSphereDyadic <- function( n, k, start="diamond", p, positive.only ) {
# Generate an approximately uniform mesh on the unit ball (in p-norm) 
# in R^n, based on repeated 2 subdivision of a starting shape, which is
# specified by start.
#
# Input:
#   n = dimension of space
#   k = number of times the starting simplex is recursively subdivided
#   start="diamond" (default) is the "diamond" formed by the unit simplex 
#     and it's rotations.  
#   start="icos3d" is the 3-dim. icosahedron with 12 vertices and 20 equal 
#     area triangular faces.
#   p = power used to compute the ell-p norm
#   positive.only=TRUE to get only the first quadrant/octant/orthant, = FALSE to get whole ball
#

n <- as.integer(n)
k <- as.integer(k)
stopifnot( n > 1, k >= 0, p > 0, is.logical(positive.only), start %in% c("diamond","icos3d" ) )

# calculate the starting simplex
if (start=="diamond") {
  # use the unit simplex and it's reflections/rotations
  a <- UnitSphereEdgewise( n=n, k=1, p=2, positive.only=positive.only )
} else {
  if (start=="icos3d") {
    if (n == 3) {
      if (!positive.only) { a <- Icosahedron( ) } 
      else { stop("start='icos3d' not allowed when positive.only=TRUE") }
    } else {
      stop('Starting shape "icos3d" is only allowed in 3-dimensions')
    }
  }
}

V <- a$V
SVI <- a$SVI    

# recursively subdivide initial simplices
if (k > 0) {
  for (i in 1:k) {
    b <- EdgeSubdivisionMulti( V, SVI, k=2, normalize=TRUE, p=p )
    V <- b$V
    SVI <- b$SVI
  }  
}

# compute all the simplices for use by plotting and integration routines
nS <- ncol(SVI)
S <- array( 0.0, dim=c(n,n,nS))
for (j in 1:nS) {
  S[,,j] <- V[SVI[,j], ]
}

mesh <- list( type="UnitSphere, dyadic", mvmesh.type=4L, n=n, m=n-1L, vps=n, S=S, V=V, SVI=SVI, 
           k=k, p=p, start=start, positive.only=positive.only ) 
class(mesh) <- "mvmesh"
return(mesh) }
#####################################################################
UnitBall <- function( n, k=1, method="dyadic", p=2, positive.only=FALSE ) {
# Generate a simplicial approximation to the unit ball consisting
# of simplices with the origin at one vertex and the vertices of 
# hyperspherical triangles on the unit sphere as the other vertices.
#
# Output is an object of class "mvmesh"
#

n <- as.integer( n )
k <- as.integer( k )
sphere <- UnitSphere( n=n, k=k, method=method, p=p, positive.only=positive.only )

origin <- rep(0.0,n)
# not sure why deparse.level is needed here, but otherwise it gives messy dimnames
newV <- rbind( sphere$V, origin, deparse.level=0 )  
nV <- nrow(newV)
nS <- dim(sphere$S)[3]
newS <- array( 0.0, dim=c(n+1,n,nS) )
for (i in 1:nS) {
  newS[,,i] <- rbind( sphere$S[,,i], origin )
}
newSVI <- rbind( sphere$SVI, rep( nV, nS ) )

mesh <- list(type=paste("UnitBall,",method),mvmesh.type=ifelse(method=="edgewise",5L,6L),
          n=n,m=n,vps=n+1L,S=newS,V=newV,SVI=newSVI,
          k=k,method=method,p=p,positive.only=positive.only)
class(mesh) <- "mvmesh"
return(mesh) }
#####################################################################
PolarSphere <- function( n, breaks=c(rep(4,n-2),8), p=2, positive.only=FALSE ) {
# Generate a simplicial approximation to the unit sphere constructed using
# the polar grid, i.e. the image of a rectangular grid on the angle space
# theta.
#
# Output is an object of class "mvmesh"
#

n <- as.integer( n )
stopifnot( n >= 2 )
# calculate rectangular mesh on angle space
rect.mesh <- RectangularMesh( a=rep(0.0,n-1), b=c(rep(pi,n-2),2*pi), breaks=breaks)

# use polar coordinate transformation to find mesh on sphere
nV <- nrow(rect.mesh$V)
V <- Polar2Rectangular( r=rep(1.0,nV), theta=rect.mesh$V )
nS <- dim(rect.mesh$S)[3]
vps <- as.integer(2^(n-1))

# construct the list of simplices by inherited grouping in rect.mesh
S <- array( 0.0, dim=c(vps,n,nS) )
for (k in 1:nS) {
  for (j in 1:vps) {
    S[j,,k] <- V[ rect.mesh$SVI[j,k],  ]
  }
}

mesh <- list(type="PolarSphere",mvmesh.type=9L,
          n=n,m=n-1L,vps=vps,S=S,V=V,SVI=rect.mesh$SVI,
          breaks=breaks,p=p,positive.only=positive.only)
class(mesh) <- "mvmesh"
return(mesh) }
#####################################################################
PolarBall <- function( n, breaks=c(rep(4,n-2),8), p=2, positive.only=FALSE ) {
# Generate a simplicial approximation to the unit ball consisting
# of simplices with the origin at one vertex and the vertices of 
# hyperspherical triangles on the unit sphere as the other vertices.
#
# Output is an object of class "mvmesh"
#

n <- as.integer( n )
sphere <- PolarSphere( n=n, breaks=breaks, p=p, positive.only=positive.only )

origin <- rep(0.0,n)
# not sure why deparse.level is needed here, but otherwise it gives messy dimnames
newV <- rbind( sphere$V, origin, deparse.level=0 )  
nV <- nrow(newV)
nS <- dim(sphere$S)[3]
vps <- as.integer(2^(n-1)+1)
newS <- array( 0.0, dim=c(vps,n,nS) )
for (k in 1:nS) {
  newS[,,k] <- rbind( sphere$S[,,k], origin )
}
newSVI <- rbind( sphere$SVI, rep( nV, nS ) )

mesh <- list(type="PolarBall",mvmesh.type=10L,
          n=n,m=n,vps=vps,S=newS,V=newV,SVI=newSVI,
          breaks=breaks,p=p,positive.only=positive.only)
class(mesh) <- "mvmesh"
return(mesh) }
#####################################################################
RectangularMesh <- function( a, b, breaks=5, silent=FALSE ) {
# Generate a mesh on the hyperrectangle [a,b]
# breaks  specifies the bin boundaries in each coordinate, one of:
#  - a vector of length d, which specifies the number of bins in each dimension. Then
#    the bin boundaries in coordinate j are given by seq(min(x[,j]),max(x[,j]),length=breaks[j])
#  - a single number m, which is equivalent to breaks=rep(m,d)
#  - a list with breaks[[j]] a vector of doubles that gives the
#    dividing points/bin boundaries for the j-th coordinate of x.  In this case, 
#    the bounds a and b are ignored.
#    In all cases, breaks is converted to a list as in the last case
# silent indicates whether or not a warning is printed if the subdivision specified by
#    breaks does not coincide with the hyperrectangle [a,b]
#
# Output is a list of class "mvmesh"
#
stopifnot( is.numeric(a), is.numeric(b), is.vector(a), is.vector(b), 
           length(a)==length(b), length(a) >= 1, all(a < b) )
n <- length(a)
orig.breaks <- breaks

# if breaks is a numeric vector, construct a list with the full set of dividing points
if ( is.numeric(breaks) & is.vector(breaks) ) {
  breaks <- as.integer(breaks)
  if( length(breaks)==1 ) { 
    breaks <- rep(breaks,n) 
  }
  stopifnot( all(breaks > 0), length(breaks)==n )
  new.breaks <- vector(mode="list",length=n)
  for (j in 1:n) {
    new.breaks[[j]] <- seq( a[j], b[j], length=breaks[j]+1)
  } 
  breaks <- new.breaks
}

# check that breaks is a valid list
if( !is.list(breaks) ) { stop("breaks must be an integer, a vector of integers, or a list") }
for (j in 1:n) {
  stopifnot( is.numeric(breaks[[j]]), is.vector(breaks[[j]]), length(breaks[[j]]) > 1, 
             all(diff(breaks[[j]])>0) )
}

if( !silent ) {
  # check to see if the specified breaks do NOT cover the region [a,b]
  aa <- rep(NA,n)
  bb <- rep(NA,n)
  for (j in 1:n) {
    tmp <- range( breaks[[j]] )
    aa[j] <- tmp[1]
    bb[j] <- tmp[2]
  }
  if ( any(a < aa) | any(b > bb) ){
    warning("specified breaks do not cover the hyperectangle [a,b]" )
  }
}

# Now construct the mesh
# number of dividing points in each coordinate
ngrid <- as.integer(sapply(breaks,length)) # as.integer( ) strips off names
# number of bins in each coordinate
nbins <- ngrid - 1L

# compute the vertices/grid points
nV <- prod(ngrid)
V <- matrix( 0.0, ncol=n, nrow=nV )
i <- rep(1L,n)
j <- 1L
while (i[1] > 0) {
  for (l in 1:n) { V[j,l] <- breaks[[l]][i[l]] }
  j <- j + 1L
  i <- NextMultiIndex(i,ngrid)
}

# compute the simplices
nS <- prod(nbins)
vps <- as.integer(2^n) # number of vertices per simplex (hyperrectangle)
S <- array( 0.0, dim=c(vps,n,nS) )
SVI <- matrix( 0L, nrow=vps, ncol=nS )
i <- rep(1L,n)
j <- 1L
while (i[1] > 0) {
  for (k in 1:vps) {
    m <- ConvertBase( k-1L, 2L, n )
    for (l in 1:n) {
      S[k,l,j] <- breaks[[l]][i[l]+m[l]]
    }
    SVI[k,j] <- MatchRow( S[k,,j], V )
  }
  i <- NextMultiIndex(i,nbins)
  j <- j + 1L
}

mesh <- list(type="RectangularMesh",mvmesh.type=7L,n=n,m=n,vps=as.integer(2^n),S=S,V=V,SVI=SVI,
     a=a,b=b,breaks=breaks)
class(mesh) <- "mvmesh"
return(mesh) }
#####################################################################
HollowTube <- function( n, k.x=1, k.circumference=2, method="dyadic", p=2 ){
# define a 'horizontal' hollow tube in n dimesions, an (n-1) dimensional 
#    surface of SolidTube( )

n <- as.integer(n)
k.x <- as.integer(k.x)
k.circumfernce <- as.integer(k.circumference)
stopifnot( n >= 2, k.x >= 1, k.circumference >= 1 )

x.grid <- seq(0,1,length=k.x+1)
vps <- 2*(n-1)
if (n==2){
  # in two dimensions, treat each line segment as a separate, unconnected 1-simplex
  nS <- 2*k.x
  SVI <- matrix( 0L, nrow=vps, ncol=nS )
  V <- cbind( rep(x.grid,2), c(rep(1,k.x+1),rep(-1,k.x+1) ) )
  for (i in 1:k.x) {
    SVI[ ,i] <- c(i,i+1)
    SVI[ ,i+k.x] <- c(i+k.x+1,i+k.x+2)    
  }
} else {
  # n > 2, use cross product of x.grid and (n-1)-dim. sphere
  sphere <- UnitSphere(n=n-1,k=k.circumference,method=method,p=p)
  nS.sphere <- dim(sphere$S)[3] # number of simplices in sphere
  nV.sphere <- nrow(sphere$V)

  V <- matrix( 0.0, nrow=(k.x+1)*nV.sphere, ncol=n )
  k <- 0
  for (i in 1:(k.x+1)) {
    for (j in 1:nV.sphere) {
      k <- k + 1
      V[k, ] <- c( x.grid[i], sphere$V[j,] )
    }
  }
  nS <- k.x * nS.sphere  
  SVI <- matrix( 0L, nrow=vps, ncol=nS )  
  k <- 0
  for (i in 1:k.x) {
    for (j in 1:nS.sphere) {
      k <- k + 1
      SVI[ ,k] <- c( (i-1)*nV.sphere + sphere$SVI[,j], rev(i*nV.sphere + sphere$SVI[,j]) ) 
    }
  }  
}

# define S array using V and SVI
S <- array( 0.0, dim=c(vps,n,nS ) )
for (i in 1:nS) {
  for (j in 1:vps) {
    S[j, ,i] <- V[SVI[j,i],]
  }
}

mesh <- list(type="HollowTube",mvmesh.type=11L,n=n,vps=vps,S=S,V=V,SVI=SVI,
    k.x=k.x,k.circumference=k.circumference)
class(mesh) <- "mvmesh"
return(mesh)}
#####################################################################
SolidTube <- function( n, k.x=1, k.circumference=2, 
    method="dyadic", p=2 ){
# define a 'horizontal' solid rod in n dimensions, an n dimensional 
# solid region inside HollowTube( )

n <- as.integer(n)
k.x <- as.integer(k.x)
k.circumfernce <- as.integer(k.circumference)
stopifnot( n >= 2, k.x >= 1, k.circumference >= 1 )

x.grid <- seq(0,1,length=k.x+1)
vps <- 2*n
if (n==2){
  # in two dimensions, treat as 2-simplex
  nS <- k.x
  SVI <- matrix( 0L, nrow=vps, ncol=nS )
  V <- cbind( rep(x.grid,2), c(rep(1,k.x+1),rep(-1,k.x+1) ) )
  for (i in 1:nS) {   
    SVI[ ,i] <- c(i,i+1,i+k.x+2,i+k.x+1)    
  }
} else {
  # n > 2, use cross product of x.grid and (n-1)-dim. ball
  ball <- UnitBall(n=n-1,k=k.circumference,method=method,p=p)
  nS.ball <- dim(ball$S)[3] # number of simplices in ball
  nV.ball <- nrow(ball$V)

  V <- matrix( 0.0, nrow=(k.x+1)*nV.ball, ncol=n )
  k <- 0
  for (i in 1:(k.x+1)) {
    for (j in 1:nV.ball) {
      k <- k + 1
      V[k, ] <- c( x.grid[i], ball$V[j,] )
    }
  }
  nS <- k.x * nS.ball  
  SVI <- matrix( 0L, nrow=vps, ncol=nS )  
  k <- 0
  for (i in 1:k.x) {
    for (j in 1:nS.ball) {
      k <- k + 1
      SVI[ ,k] <- c( (i-1)*nV.ball + ball$SVI[,j], rev(i*nV.ball + ball$SVI[,j]) ) 
    }
  }   
}

# define S array using V and SVI
S <- array( 0.0, dim=c(vps,n,nS ) )
for (i in 1:nS) {
  for (j in 1:vps) {
    S[j, ,i] <- V[SVI[j,i],]
  }
}


mesh <- list(type="SolidTube",mvmesh.type=12L,n=n,vps=n,S=S,V=V,SVI=SVI,
    k.x=k.x,k.circumference=k.circumference)
class(mesh) <- "mvmesh"
return(mesh)}
#####################################################################
NextMultiIndex <- function( i, n ) {
# compute the next multi-index beyond (i[1],...,i[k]), with each i[j] in the range 1:n[j]
# it is assumed length(i)=length(n).  Essentially, we count with different bases for
# each position

cur <- 1L
k <- length(i)
repeat {
  i[cur] <- i[cur]+1L
  if (i[cur] <= n[cur]) { break }
  if (cur == k) {
    i[1] <- 0
    break
  }
  i[cur] <- 1L
  cur <- cur + 1L
}
return(i) }
#####################################################################
EdgeSubdivision <- function( n, k ){
#
#  This is the core routine for subdivision, which 
#  implements the algorithm in Edgewise subdivision of a simplex, 
#  by Edelsbrunner and Grayson, Discrete Comput. Geom., Vol 24, 707-719 (2000).
#
#  This function computes the "color scheme" for a  equal volume 
#  subdivision of an n-dimensional simplex.  This color scheme is a 
#  "base k coding" of the vertices of the subsimplices.  The coding is independent
#  of a simplex; it is a listing of the "base k barycentric coordinates"
#  that can then be used on any simplex, including one that sits inside a higher 
#  dimensional space.   EdgeSubdivision is based on matlab code by 
#  Goncalves, Palhares, Takahashi, and Mesquita, Algorithm 860: SimpleS -- 
#  an extension of Freudenthal's simplex subdivision, ACM Trans. Math. Softw., 
#  32, 609-621 (2006).
#
# Input:
#   n = dimension of the simplex
#   k = number of subdivisions of each edge
#
# Output: 
#   T = an array of size  k x (n+1) x k^n,  where each T[,,i] is
#       a k x (n+1) "color scheme" matrix, which determines a sub-simplex
#
#  This routines is "coordinate free", it works with a color scheme, that 
#  describes a way of  subdividing a simplex. So the output array T
#  is independent of point coordinates, it can be
#  used to find a uniform mesh for any (every) n dimensional simplex 
#  inside R^(n+p), p >= 0.  Use function SimplexCoord( ) to get 
#  the coordinates of a subsimplex.

n <- as.integer(n)
k <- as.integer(k)
stopifnot( length(n)==1, length(k)==1 )

nS <- k^n
T <- array( 0L, dim=c(k,n+1,nS) )
CS <- matrix( 0L, nrow=k, ncol=n+1)  

for (l in 1L:nS) {
    x <- ConvertBase( l-1L, k, n )
    cor <- 0L
    for (i in 1L:k){
        CS[i,1] <- cor
        for (j in 1:n) {
            if ( x[n+1-j] == (i-1L) ) { cor <- cor + 1L }
            CS[i,j+1] <- cor;
        }
    }
    T[,,l] <- CS
}
return(T)  }
#####################################################################
ConvertBase <- function( m, b, n ){
# convert the positive integer x to base b, using n digits:
#   m  = y[1] + b*y[2] + b^2*y[3] + ... + b^(n-1)*y[n]

y <- rep(0L,n)  
for(i in 1:n){
  d <- floor(m/b)
  y[i] <- m-d*b
  m <- d
}
return(y) }        
#####################################################################
SimplexCoord <- function( S, color ){
# compute the coordinates of one m-dim. sub-simplex of S with color scheme chi
#
# Input:
#   S    vps x n matrix, with the rows giving the vertices of a simplex in R^n 
#   color  k x vps color scheme matrix
#
# Output: 
#   P   vps x n matrix, with rows giving the coordinates of sub-simplex of S 

n <-ncol(S)
vps <- nrow(S)  # m+1
P <- matrix(0.0,nrow=vps,ncol=n)
for (j in 1:vps) {
  P[j,] <- PointCoord(S,color[,j])
}
return(P) }
#####################################################################
PointCoord <- function( S, color ){
# compute the coordinates of the point in simplex S with color 
# scheme in 'color'.  color/k gives the barycentric coordinates of 
# the point inside S

n <- ncol(S)
k <- length(color)
P <- rep( 0.0, n )
for (i in 1:k) { P  <- P + S[color[i]+1, ]   }
return(P/k) }
#####################################################################
SVIFromColor <- function( S, T ) {
#
# Find the SVI of a single simplex S determined by the color schemes in array T
#
# Input
#    S is a vps x n matrix, with S[1,],...,S[vps,] giving the 
#      points in R^n that are the vertices of the simplex
#    T is a k x n x nS array with T[,,j] being the color scheme for
#      the j-th subsimplex of S
# 
# Output is a list with 2 fields:   
#    V an nV x n matrix; with rows V[1,],..,V[nV,] giving the 
#         vertices/coordinates in R^n
#    SVI an integer matrix of size n x nS, with column SVI[,1],...,SVI[,L]
#         giving the indices for rows of V that make up each subsimplex

n <- ncol(S)
k <- dim(T)[1]
vps <- dim(T)[2]  
nS <- dim(T)[3]

# find the vertices, eliminating duplicates
# put all columns of all color schemes into one large matrix
a <- matrix( 0L, nrow=k, ncol=nS*vps )
m <- 0
for (i in 1:nS) {
  for (j in 1:vps) {
    m <- m + 1
    a[,m] <- T[,j,i]
  }
}

# eliminate duplicate columns 
b <- unique(a,MARGIN=2)  

# compute vertices of each (now unique) color scheme row
nV <- ncol(b)
V <- matrix(0.0,nrow=nV,ncol=n)
for (i in 1:nV) {
  V[i,] <- PointCoord(S,b[,i])
}

# compute pointers to vertices
SVI <- matrix( 0L, nrow=vps, ncol=nS )
for (i in 1:vps) {
  for (j in 1:nS) {
    for (m in 1:nV) {
      if( identical(b[,m],T[,i,j] )) { SVI[i,j] <- m }
    }
  }
}  

return( list( V=V, SVI=SVI ) ) }
#######################################################################
NumVertices <- function( n, k, single=TRUE ) {
# compute v[k,n] = number of vertices in k-th order edge subdivision of 
#    a (n-1)-dimensional simplex in R^n
# if single=TRUE, return a single value v[k,n]
# otherwise, return the whole matrix v[1:k,1:n]

# use a recursion formula
  v <- matrix(as.integer(0),k,n)
  v[1,] <- 1:n
  v[,1] <- as.integer(1)
  for (j in 2:n) {
    for (i in 2:k) {
      v[i,j] <- v[i-1,j] + v[i,j-1]
    }
  }
  if (single) return(v[k,n])
  else return(v)
}
#####################################################################
MatchRow <- function( v, table, first=1, last=nrow(table) ) {
# Search through rows first, first+1,...,last of table to see 
# if v is present. If so, return row numbers of ALL matches

j <- integer(0)
for (i in first:last) {
  if (identical(v,table[i,]) ) { j <- c(j,i) }
}

return( j ) }
#####################################################################
EdgeSubdivisionMulti <- function( V, SVI, k, normalize=FALSE, p=2 ) {
# Do an edge k-subdivision of multiple simplices, guaranteeing
# that the result has no repeated vertices.  This function was written
# for UnitSphereDyadic; it is not general.
#
# Input:
#   V = nV x n matrix of vertices, rows V[1,],...,V[nV,] are vertices in R^n
#   SVI = vps x nS matrix, with columns SVI[,1],...,SVI[,nSV] specifying
#       a vps-dim. simplex, i.e. indices of the vertices in V
#   normalize=TRUE means the vectors in V are normalized to have unit length
#   p = the norm to use if normalizing, e.g. p=2 is Euclidean norm
#
# Output:
#   V = nVnew x n matrix of vertices, guaranteed to have no repeats
#   SVI = vps x nSnew matrix giving the grouping
#

k <- as.integer( k )
stopifnot( is.matrix(V), is.matrix(SVI), k > 1, is.logical(normalize), p > 0.0 )
n <- ncol(V); nV <- nrow(V); vps <- nrow(SVI); nS <- ncol(SVI)

# compute first subdivision
T <- EdgeSubdivision( n=vps-1, k=k )
new <- SVIFromColor( V[ SVI[,1], ], T )
M <- nrow(new$V); L <- ncol(new$SVI)

# using sizes from first subdivision, allocate matrices for subdivisions
Vnew <- matrix(0.0,ncol=n,nrow=nS*M)
SVInew <- matrix( 0L, nrow=vps, ncol=nS*L)

# loop through simplices, subdividing each
Vend <- 0L
SVIend <- 0L
for (i in 1:nS) { 
  new <- SVIFromColor( V[SVI[,i], ], T )
  Vnew[ (Vend+1):(Vend+M), ] <- new$V
  SVInew[ ,(SVIend+1):(SVIend+L)] <- new$SVI + Vend
  Vend <- Vend+M    
  SVIend <- SVIend + L
}


# eliminate duplicate vertices. Duplicates happen when some of original
# simplices share edges
count <- 1
for (i in 2:Vend) {
  l <- which(SVInew==i)  
  j <- MatchRow( Vnew[i,], Vnew, 1, count )
  if (length(j)==0) { 
    # Vnew[,i] has not occurred yet in the matrix, insert it 
    count <- count + 1
    Vnew[count,] <- Vnew[i,]  
    SVInew[l] <- count
  } else { 
    SVInew[l] <- j[1]       
  }
}

# If requested, normalize to be on unit sphere in l^p norm
if (normalize) {
  r <- LpNorm(Vnew[1:count,],p)
  for (i in 1:count) {
    Vnew[i,] <- Vnew[i,]/r[i]
  }
}  

mesh <- list( V=Vnew[1:count,], SVI=SVInew) 
class(mesh) <- "mvmesh"
return(mesh)}
#####################################################################
Icosahedron <- function( ) {
# return the vertices and list of simplices for an icosahedron in 
#   R^3 with vertices on the unit sphere.  12 vertices, 20 faces

# specify the coordinates of the icosahedron, then normalize to have unit length
p <- (sqrt(5)-1)/2  # phi="little" golden ratio
V <- matrix( c(0,p,1,  0,p,-1,  0,-p,1,  0,-p,-1,  p,1,0,  p,-1,0, 
               -p,1,0, -p,-1,0, 1,0,p,  -1,0,p,   1,0,-p, -1,0,-p), nrow=12, ncol=3, byrow=TRUE)
c1 <- 1.0/sqrt( sum( V[1,]^2 ) )
V <- c1*V

# specify the grouping matrix
SVI <- matrix( as.integer(c( 1,3,9,  1,9,5,  1,5,7,  1,7,10,  1,10,3,  4,12,2, 
    4,2,11,  4,11,6,  4,6,8,  4,8,12,  9,3,6,  5,9,11,  7,5,2,  10,7,12,  
    3,10,8,  2,12,7,  11,2,5,  6,11,9,  8,6,3,  12,8,10 ) ), nrow= 3 , ncol= 20 )
S <- array( 0.0, dim=c(3,3,20))
for (k in 1:20) {
  S[ ,,k] <- V[ SVI[,k],  ]
}
mesh <- list( type="Icosahedron", mvmesh.type=8L, n=3L,k=NA, originalS=S, S=S, V=V, SVI=SVI )
class(mesh) <- "mvmesh"
return( mesh ) }
#####################################################################
LpNorm <- function( x, p ){ 
# compute the l^p norm of the rows of x: |x|_p = ( sum(abs(x)^p) )^(1/p)
# When p=Inf, it correctly returns max( abs(x) ).

if (is.vector(x)) { x <- matrix(x,nrow=1) }
stopifnot( is.matrix(x), nrow(x) > 0, ncol(x) > 0, p > 0 )

if (is.finite(p)) {
  r <- apply( abs(x)^p, MARGIN=1, FUN=sum )^(1/p)
} else {
  r <- apply( abs(x), MARGIN=1, FUN=max )
}
return(r) }
###################################################################
AffineTransform <- function( mesh, A, shift ) {
# form new mesh by tranforming each vertex v to A %*% v + shift

stopifnot( is.list(mesh), class(mesh)=="mvmesh", is.matrix(A), is.numeric(A), is.numeric(shift),
   is.vector(shift), nrow(A)==ncol(A), nrow(A)==length(shift) )

# copy to new mesh and then modify the appropriate fields   
new <- mesh
new$type <- paste( mesh$type, "+affine", sep="" )
nV <- nrow( new$V )
nS <- ncol( new$SVI )
vps <- new$vps

# transform vertices
affine.func <- function( v, A, b ) { A %*% v + b }
new$V <- t( apply( mesh$V, MARGIN=1, affine.func, A=A, b=shift ) )

# copy coordinates of new simplices
for (k in 1:nS) {
  for (j in 1:vps) {
    new$S[j,,k] <- new$V[ new$SVI[j,k], ]
  }
}
return( new ) }
#######################################################################
Rotate2D <- function( theta ) {
# 2D rotation matrix for use by AffineTransform

stopifnot( is.numeric(theta), length(theta)==1 )

return( matrix( c( cos(theta),sin(theta),-sin(theta),cos(theta) ),
                      nrow=2,ncol=2) )}
#######################################################################
Rotate3D <- function( theta ) {
# 3D rotation matrix for use by AffineTransform

stopifnot( is.numeric(theta), length(theta)==3 )

# compute the rotations around each axis, then return their product
Rx <- matrix( c(1,0,0,  0,cos(theta[1]),sin(theta[1]),  
              0,-sin(theta[1]),cos(theta[1]) ), nrow=3, ncol=3)
Ry <- matrix( c(cos(theta[2]),0,-sin(theta[2]),  0,1,0,  
              sin(theta[2]),0,cos(theta[2]) ), nrow=3, ncol=3)
Rz <- matrix( c(cos(theta[3]),sin(theta[3]),0, -sin(theta[3]),cos(theta[3]),0,
              0,0,1 ), nrow=3, ncol=3)
return( Rx %*% Ry %*% Rz )}
#######################################################################
V2Hrep <- function( S ) {
# Convert a (vertex) V-representation for a list of simplices in array S into a 
# list H of the (half-space) H-representation of the simplices. It
# is implicit that all the simplices in S have the same number of
# vertices, and we assume that their H-representations are all of 
# the same size; an error will occur if they are not.
#
# S is an (vps x n x nS) array of simplices, S[i,,k] is the i-the vertex of the k-th simplex

if (is.matrix(S) ) { S <- array( S, dim=c(nrow(S),ncol(S),1) ) }
nS <- dim(S)[3]

dim.Hrep <- c(0L,0L)
for (k in 1:nS) {
  Vk <- makeV( S[,,k] )
  tmp <- rcdd::scdd( Vk, representation="V" )$output  
  if (k==1) { 
    dim.Hrep <- dim(tmp)
    H <- array( 0.0, dim=c(dim.Hrep,nS) ) 
  }
  if ( any(dim(tmp) != dim.Hrep)) stop( "simplices have different size H-representations" )
  H[,,k] <- tmp
}
return(H) }
#######################################################################
H2Vrep <- function( H ) {
# convert an H (half-space) representation for a polytope into a V (vertex)
# representation 
# H = nC x (n+2) x nS array of H-representations.  H[,,k] is the H-rep for
# the k-th simplex.  See the documentation for rcdd for the format of H.
# It is assumed that all H-reps have the same number of constraints, and
# that the resulting V-reps have the same number of vertices per simplex (vps)


if (is.matrix(H) ) { H <- array( H, dim=c(nrow(H),ncol(H),1) ) }
dimH <- dim(H)
stopifnot( is.array(H), length(dimH)==3 )
nS <- dimH[3]
dim.Vrep <- c(0L,0L)
for (k in 1:nS) {
  tmp <- rcdd::scdd( H[,,k], representation="H" )$output[ ,-c(1,2),drop=FALSE]
  if (k==1) {
    dim.Vrep <- dim(tmp)
    S <- array( 0.0, dim=c(dim(tmp),nS)) 
  }
  if ( any(dim(tmp) != dim.Vrep)) stop( "simplices have different size V-representations" )
  S[,,k] <- tmp
}

return(S) }
#######################################################################
SatisfyHrep <- function( x, Hsingle ) {
# determine if the vectors in x[,1],x[,2],... satisfy the constraints
# in Hsingle, an H-representation of a single simplex
# The result is a vector of integers that describe which columns satisfy 
# the constraints.

stopifnot( is.matrix(Hsingle) )
l <- Hsingle[ ,1]
b <- Hsingle[ ,2]
A <- - Hsingle[ , -(1:2)]
if (is.vector(x)) { x <- matrix( x, nrow=1 ) }
stopifnot( is.matrix(x), ncol(x)+2==ncol(Hsingle) )
n <- ncol(x)
nx <- nrow(x)

satisfy <- integer(0)
for (k in 1:nx) {
  axmb <- A %*% x[k,] - b
  if (all(axmb <= 0.0) & all( l*axmb == 0 ) ) { satisfy <- c(satisfy,k) }
}
return(satisfy)}
########################################################################
Polar2Rectangular <- function( r, theta ) {
# Convert polar coordinates to rectangular coordinates in n-dimensions.
# If r is a scalar, then the polar point (r,theta[1:(n-1)]) is converted to
#     rectangular coordinates x[1:n].
# If r is a vector of length m, then theta should be a matrix of dimensions m x (n-1),
#     and the result is a matrix x[1:m,1:n], with rows of x being the rectangular
#     coordinates of points (r[j],theta[,j]).
# The result is always a matrix x of size (m x n).
#
m <- length(r)
if (!is.matrix(theta)) { theta <- as.matrix(theta,nrow=1) }
stopifnot( m == nrow(theta))
n <- ncol(theta) + 1L
x <- matrix(0.0,nrow=m,ncol=n)
for (j in 1:m) {
  col.cos <- cos(theta[j,])
  col.sin <- sin(theta[j,])
  s <- c( col.cos[1], rep(col.sin[1],n-1) )
  if (n > 2) {
    for (k in 2:(n-1)) {
      s[k] <- s[k]*col.cos[k]
      s[(k+1):n] <- s[(k+1):n]*col.sin[k]
    }
  }
  x[j,] <- r[j]*s
}
return(x) }
########################################################################
Rectangular2Polar <- function( x ) {
# Convert from rectangular to polar coordinates in n dimensions.
# If x[1:n] is a vector in rectangular coordinates, convert to
#     spherical coordinates (r,theta[1:(n-1)])
# If x[1:m,1:n] is a matrix, convert each row x[i, ] to spherical coordinates
#   given by r[i] and theta[,i]
# Output is a list with entries
#   r[1:m] a vector of lengths
#   theta[1:m,1:(n-1)] matrix with theta[j,] giving angles for point x[j, ]

if(!is.matrix(x)) { x <- as.matrix(x,nrow=1) }
n <- ncol(x)
m <- nrow(x)
r <- rep(0.0,m)
theta <- matrix(0.0,ncol=n-1,nrow=m)
for (j in 1:m) {
  rsq <- x[j,]^2
  cum.rsq <- cumsum(rev(rsq))
  r[j] <- sqrt( cum.rsq[n] )
  if (r[j] > 0.0) {
    if (n>2) {
      for (k in 1:(n-2)) {
        theta[j,k] <- atan2( sqrt(cum.rsq[n-k]), x[j,k] )
      }
    }
    theta[j,n-1] <- 2*atan2( x[j,n], x[j,n-1]+sqrt(cum.rsq[2] ) )
  }
}
return(list(r=r,theta=theta))}
##################################################################################
IntersectMultipleSimplicesV <- function( S1, S2 ) {
# Find the pairwise intersection of two lists of simplices, both given in the V-representation.
# only the intersections with positive volume are returned
# S1 is an (n x m1 x nS1) array, with S1[,,k] being the V-representation of the k-th 
#    simplex in the first set of simplices
# S2 is an (n x m2 x nS2) array, with S2[,,k] being the V-representation of the k-th 
#    simplex in the second set of simplices
# Note that all simplices in S1 must have the same number of vertices, and all
#   simplices in S2 must have the same number of vertices.  The function Intersect2SimplicesH( )
#   handles the intersection of two arbitrary simplices, but for convenience here we
#   require unform sizes within S1 and within S2.
#   
# return a list with 3 fields:
#   S  an (n x (n+1) x nS) array, list of resulting simplices
#   index1[1:nS] integer array 
#   index2[1:nS] integer array
#   S[,,k] is in the intersection of S1[,,index1[k]] and S2[,,index2[k]]
#  Since the intersection may be an n-simplex with more than (n+1) vertices,
#  each intersection is 'triangulated' to be decomposed into simplices with exactly (n+1) vertices

H1 <- V2Hrep( S1 )
H2 <- V2Hrep( S2 )
return( IntersectMultipleSimplicesH( H1, H2 ) ) }
##################################################################################
IntersectMultipleSimplicesH <- function( H1, H2 ) {
# Similar to IntersectSimplicesV, but now the simplices are given in
# the H-representation

stopifnot( ncol(H1)==ncol(H2) )
n <- ncol(H1)-2
count <- 0
S <- array( 0.0, dim=c(n+1,n,0) )
index1 <- integer(0)
index2 <- index1
nS1 <- dim(H1)[3]
nS2 <- dim(H2)[3]
for (i2 in 1:nS2) {
  for (i1 in 1:nS1) {
    a <- Intersect2SimplicesH( H1[,,i1], H2[,,i2], tesselate=TRUE )
    if(!is.null(a$S)) { 
      m <- dim(a$S)[3]
      index1 <- c(index1,rep(i1,m))
      index2 <- c(index2,rep(i2,m))
      S <- abind( S, a$S, force.array=TRUE )      
    }
  }
}
return( list( S=S, index1=index1, index2=index2 ) )}
##################################################################################
Intersect2SimplicesH <- function( H1, H2, tesselate=FALSE ) {
# intersect two arbitrary simplices given by their H-represenations H1 and H2
# if tesselate=TRUE, the resulting set is tesselated so that we get 
#
# Return a list with fields:
#   H = H-representation of the intersection simplex
#   V = vertices (V-rep.) of the intersection simplex (NULL if intersection has zero volume)
#   S = tesselation of the intersection simplex (if tesselate=TRUE, NULL if intersection has zero volume)
cat("Intersect2...\n")
n <- ncol(H1)-2
H <- rcdd::redundant( rbind( H1, H2 ), representation="H")$output
V <- NULL
S <- NULL
if( nrow(H) > n) {
  V <- H2Vrep( H )[,,1] 
  if (tesselate & (nrow(V) > n)) {
cat("V="); print(V)
    del.tess <- geometry::delaunayn( V, options="Qz" )
    S <- array( 0.0, dim=c(n+1,n,nrow(del.tess)) )
    for (k in 1:nrow(del.tess)) {
      b <- del.tess[k,]
      cur.tess <- V[b,]
      S[,,k] <- cur.tess
    }
  }
}
return( list(H=H, V=V, S=S) ) }
##################################################################################
Lift2UnitSimplex <- function( S ){
# In several places a simplex S on the unit simplex in R^n is projected onto
#   S[,-n] in R^(n-1).  This function is the inverse of that projection:
#   it lifts that projection back to the unit simplex.

n <- nrow(S)
if (is.matrix(S)) { S <- array(S,dim=c(nrow(S),ncol(S),1) )}
stopifnot( is.array(S), length(dim(S))==3 )

nS <- dim(S)[3]
S2 <- array( 0.0, dim=c(n,n,nS) ) 
for (k in 1:nS) { 
  tmp <- 1-rowSums( S[,,k,drop=FALSE] )
  S2[,,k] <- cbind(S[,,k],tmp)
}
return(S2) }
##################################################################################
# replace with aperm( B, c(2,1,3) )
#TransposeMultipleMatrices <- function( B ) {
#  For the (n x m x p) array B, compute the (m x n x p) array C
#  with C[,,i] = t(B[,,i])

#if( is.matrix(B) ) { return( t(B) ) }
#dimB <- dim(B)
#C <- array( 0.0, dim=c(dimB[2],dimB[1],dimB[3]) )
#for (i in 1:dimB[3]) {  C[,,i] <- t(B[,,i]) }
#return(C) }
########################################################################
mvmeshFromSVI <- function( V, SVI ) {
# construct an mvmesh object from a tessellation: 
#  V = nV x n matrix of vertices; V[i,] is the i-th vertex
#  SVI = vps x nS matrix of integer indices SVI giving the grouping;
#        SVI[,j] gives the indices of vertices that make up simplex j
# Note: this assumes that the faces are all of dimension (vps-1)

vps <- nrow(SVI)
n <- ncol(V)
nS <- ncol(SVI)
S <- array( 0.0, dim=c(vps,n,nS) )
for (i in 1:nS) {
  for (j in 1:vps) {
    S[j,,i] <- V[SVI[j,i],]
  }
}

# construct a mvmesh object using computed quantities
a <- list(type="mvmeshFromSVI",mvmesh.type=-1,n=n,m=vps-1,vps=vps,S=S,V=V,SVI=SVI)
class(a) <- "mvmesh"
return(a) }
########################################################################
mvmeshFromSimplices <- function( S ) {
# construct an mvmesh object from an array of simplices in S

# if a single simplex, convert to an array
if (is.matrix(S)) { S <- array(S,dim=c(dim(S),1) ) }
stopifnot( is.array(S), length(dim(S))==3 )

# find unique vertices
dimS <- dim(S)
vps <- dimS[1]
n <- dimS[2]
nS <- dimS[3]

# get unique rows of S
V <- uniqueRowsFromDoubleArray( S )

# construct SVI matrix from V and S
SVI <- matrix( 0L, nrow=vps, ncol=nS)
for (i in 1:nrow(V)) {
  for (j in 1:nS) {
    k <- MatchRow( V[i,], S[,,j])
    if (length(k) > 0) {
      SVI[k,j] <- i
    } 
  }
}

k <- which(SVI==0)
if( length(k) > 0 ) warning(paste(length(k),"unmatched rows in function mvmeshFromSimplices"))

# construct a mvmesh object using computed quantities
a <- list(type="mvmeshFromSimplices",mvmesh.type=-1,n=n,m=vps-1,vps=vps,S=S,V=V,SVI=SVI)
class(a) <- "mvmesh"
return(a)}
########################################################################
uniqueRowsFromDoubleArray <- function( A ) {
# find unique rows of an array A of type double using the MatchRow function,
# which is based on built-in function identical( )

# treat a matrix as an array with 3rd dimension=1
if( is.matrix(A) ) { A <- array( A, dim=c(dim(A),1) ) }
stopifnot( is.array(A), length(dim(A))==3, all(dim(A)>0) )

nA <- dim(A)[3]
B <- matrix( 0, ncol=ncol(A), nrow=nrow(A)*nA )
B[1,] <- A[1,,1]
count <- 1
for (i in 1:nA) {
  for (j in 1:nrow(A)) {
    match <- MatchRow( A[j,,i], B, first=1, last=count )
    if (length(match)==0) { 
      count <- count + 1 
      B[count,] <- A[j,,i]
    }
  }
}
return(B[1:count, ]) }
########################################################################
mvmeshFromVertices <- function( V ) {
# construct an mvmesh object from a matrix of vertices; the constructed
# simplices are ALWAYS of dimension n=ncol(V) because delaunayn will return n-dim. simplices.
# So this will NOT work with simplices of dimension m < n, i.e. spheres or 
# lower dim. simplices.  Even when the simplices are of dimension n, you may get a 
# different tessellation than you were expecting: vertices can be grouped in multiple ways.

# tesselate
SVI <- t( delaunayn( V ) )

# define sizes
n <- ncol(V)
nV <- nrow(V)
vps <- nrow(SVI)

# construct list of simplices
S <- array( 0.0, dim=c(vps,n,ncol(SVI)))
for (i in 1:ncol(SVI)) {
  for (j in 1:n) {
    S[j,,i] <- V[SVI[j,i],]
  }
}

# construct a mvmesh object using computed quantities
a <- list(type="mvmeshFromVertices",mvmesh.type=-1,n=n,m=n,vps=vps,S=S,V=V,SVI=SVI)
class(a) <- "mvmesh"
return(a)}
######################################################################
rtessellation <- function( n, S, weights=rep(1,dim(S)[3]) ) {
# simulate n random vectors from a tessellation
# S is a (vps x d x nS) array of simplices
# weights is a vector of weights

# special case of one simplex
if( is.matrix(S) ) { S <- array(S,dim=c(dim(S),1)) }

stopifnot( n > 0, is.array(S), length(dim(S))==3, dim(S)[3]==length(weights))
dimS <- dim(S)
vps <- dimS[1]; d <- dimS[2];  nS <- dimS[3]

# generate u uniform on the unit simplex
u <- matrix( rexp( n*vps),nrow=vps,ncol=n )
u.sum <- colSums(u)
for (i in 1:n) { u[,i] <- u[,i]/u.sum[i] }

x <- matrix( 0.0, nrow=n, ncol=d )
which.splx <- sample.int( n=nS, size=n, replace=TRUE, prob=weights )
for (i in 1:n) {
  j <- which.splx[i]
  x[i,] <- u[,i] %*% S[,,j]
}
return(x) }
######################################################################
mvmeshCombine <- function( mesh1, mesh2 ) {
# combine 2 meshes, retain type of mesh1

stopifnot( mesh1$n == mesh2$n, mesh1$m==mesh2$m, mesh1$vps==mesh2$vps )
new <- mesh1
new$V <- rbind(mesh1$V,mesh2$V)
new$SVI <- cbind(mesh1$SVI,mesh2$SVI)
new$S <- abind( mesh1$S, mesh2$S )
return(new) }
######################################################################
