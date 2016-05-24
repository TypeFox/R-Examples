# Multivariate histograms.  John Nolan, jpnolan@american.edu March 2015
#
# A variety of different histotgrams are given:
#   histDirectional for directional histograms
#   histRectangular for rectangular grid histograms
#   histSimplex for simplicial histograms
#
# The histograms work in any dimension d, though the plots may be difficult to
# interpret when dimension is greater than 3. 
#
# Input values common to all histograms:
#     x is the data matrix, an (n x d) matrix with row x[j,]=j-th data vector
#        NOTE: in these functions, d=dimension of data, n=sample size. In the core mvmesh functions,
#        n is the dimension.
#     plot.type determines how the histogram is ploted.
#     freq = TRUE for a frequency histogram, FALSE for a relative frequency (probability) histogram
#     report specifies the action to take with a tie, e.g. a data point lying on 
#        a boundary that is shared by 2 or more regions.  Possible values:
#            "all" print a message for every tie and reject (=data point not in any cone)
#            "summary" print a summary message only
#            "none" don't report ties
#      ...  optional parameters like col="blue", "lwd=3", etc. for plot    
#
# Each of the types of histograms has a corresponding tally function, which
# counts how may data points lie in each region.  The number of regions and
# type of region varies with the histogram type.
#
# Output is a plot of the histogram (unless plot.type="none"), with a variety of plot types, 
# and a list with multiple fields.  In all cases, the list includes fields:
#    counts[1:m] = the number of data points in each region
#    nrejects = number of data points not in any region
#    nties = number of points on the boundary of more than one region
# Other output fields vary with type of histogram.
#
# At the lowest level, each bin is specifed by a simplex specified in the H-representation
# (Half-space representation) by a matrix H.  We use the method NOT recommended by Fukuda 
# on his FAQ in Polyhedral Computation located at  http://www.inf.ethz.ch/personal/fukudak/polyfaq/node22.html
# We do this for efficiency reasons: here the number of constraints in H is
# not large, and we will generally have n=ncol(x)=sample size much larger, so we want an efficient
# way of testing whether many points lie in a set of cones.
# 
###############################################################################
histDirectional <- function( x, k, p=2, plot.type="default", freq=TRUE,
                              positive.only=FALSE, report="summary", label.octants=TRUE, 
                              normalize.by.area=FALSE, ... ) {
# plot a directional histogram for the d-dimensional data in x
# Written by John Nolan, 25 Oct. 2014
# Function specific input:
#    k is the number of subdivisions along each edge of unit simplex
#    p is the power in the l^p norm:  |x| = (sum(abs(x)^p))^(1/p)
#    plot.type specifies the type of plot to do
#    freq TRUE for a frequency histogram, FALSE for a relatuve frequency histogram
#    positive.only=TRUE for all data is positive, =FALSE for all directions
#
# Output: a list returned invisibly containing the tally and simplex info
#  

# input checking, more checking is done in directional.tally
k <- as.integer( k )
allowed.types <- c("default", "index","radial","none","orthogonal","grayscale")
allowed.report <- c("all","summary","none")
stopifnot( k > 0, is.matrix(x), nrow(x) > 1, ncol(x) > 0,  is.logical(positive.only), is.numeric(p), length(p)==1,
   plot.type %in% allowed.types, report %in% allowed.report  ) 
d <- ncol(x)
if (plot.type == "default") {
  plot.type <- ifelse( d <= 3, "radial","index")
} 
if ( (plot.type != "index") & (d > 3) ) { 
  stop( 'when dimension of data > 3, plot.type must be "index".' )
}
if ( (plot.type == "grayscale") & (d != 3) ) { 
  stop( 'plot.type "grayscale" is only allowed when the dimension of data = 3' )
}
if ( plot.type == "orthogonal" & ( (d != 3) | !positive.only | p != 1 ) ) { 
  stop( 'plot.type "orthogonal" is only allowed when the dimension of data = 3' )
}

# find the directional simplices/generators, H-rep for the cones and then
# tally the number of data points in each direction
sphere <- UnitSphere( n=d, k=k, method="dyadic", p=p, positive.only=positive.only )
S <- sphere$S
H <- HrepCones( S )
tally <- TallyHrep( x, H, report )
if (freq) {
  value <- tally$counts
  ylab.string <- "frequency"
} else {
  value <- tally$rel.freq
  ylab.string <- "relative frequency"  
}

if (normalize.by.area) {
  for (i in 1:length(value)) {
    value[i] <- tally$counts[i]/SimplexSurfaceArea( t(S[,,i]) )
  }
  ylab.string <- "frequency/area"
}

# plot the results
if ( plot.type != "none" ) {
  ncones <- length(value)
  if (plot.type == "index") {
    val.max <- max(value)
    cone <- 1:ncones
    if (positive.only) {
      n.octants <- 1
    } else {
      n.octants <- 2^d    
    }
    octant.width <- ncones/n.octants
    plot(cone,value,type='h',ylab=ylab.string,ylim=c(0,val.max), ...)
    if (label.octants) {
      plus.minus <- c("+","-")
      for (i in 1:n.octants) {
        y <- ConvertBase(i-1,2,d)
        signs <- (-1)^y
        label <- paste( plus.minus[y+1],collapse="")
        x.center <- (i-1/2)*octant.width
        text( x.center,val.max, label )
        abline(v=i*octant.width+1/2,col='red',lty=2)
      }
    }
  } else {
    scale <- max(value)  
    if (plot.type=="radial") {  
      if (d==2) {
        # set plot window
        plot( 2*range(S[,1,]), 2*range(S[,2,]), type='n',xlab="",ylab=ylab.string) 
        for (j in 1:ncones) {
          DrawSimplex2d( S[,,j], label=j, show.labels=FALSE, mvmesh.type=sphere$mvmesh.type, ... )
          mid <- colMeans( S[,,j] )
          tmp <- 1+value[j]/scale
          lines( c(1,tmp)*mid[1], c(1,tmp)*mid[2], ... )
        }
      } else {  # d==3 
        open3d()
        for (j in 1:ncones) {
          DrawSimplex3d( S[,,j], label=j, show.labels=FALSE, mvmesh.type=sphere$mvmesh.type, ... )
          mid <- colMeans( S[,,j] )
          tmp <- c(1,1+value[j]/scale)
          lines3d( tmp*mid[1], tmp*mid[2], tmp*mid[3], ... )      
        }  
      }
    }
    if (plot.type=="orthogonal") { # must have d=3
      open3d()
      # add axes and vertex labels
      lines3d( c(0,1.5), c(0,0), c(0,0) )
      lines3d( c(0,0), c(0,1.5), c(0,0) )
      lines3d( c(0,0), c(0,0), c(0,1.5) )
      text3d( 1.2,0,0, "(1,0,0)")
      text3d( 0,1.2,0, "(0,1,0)")
      text3d( 0,0,1.2, "(0,0,1)")    
      rgl.viewpoint( theta=30, phi=30)  
      
      # plot triangles and spike   
      for (j in 1:ncones) {
        a <- S[,,j]
        DrawSimplex3d( a, show.labels=FALSE, mvmesh.type=sphere$mvmesh.type, ... )
        mid <- colMeans( a )
        tmp <- value[j]/scale  
        lines3d( c(mid[1],mid[1]+tmp), c(mid[2],mid[2]+tmp), c(mid[3],mid[3]+tmp), ... )
      }   
    }
    if (plot.type=="grayscale") { # d must be 3
      open3d()
      # add axes and vertex labels
      lines3d( c(0,1.5), c(0,0), c(0,0) )
      lines3d( c(0,0), c(0,1.5), c(0,0) )
      lines3d( c(0,0), c(0,0), c(0,1.5) )
      text3d( 1.2,0,0, texts="(1,0,0)")
      text3d( 0,1.2,0, texts="(0,1,0)")
      text3d( 0,0,1.2, texts="(0,0,1)")    
      rgl.viewpoint( theta=30, phi=30)   
      
      # next plot all the triangles, shading proportional to count in that direction    
      for (j in 1:ncones) {
        a <- S[,,j]
        DrawSimplex3d( a, show.labels=FALSE, mvmesh.type=sphere$mvmesh.type, ... )
        triangles3d( a, col=gray(1-value[j]/scale), specular="black"  )
      }       
    
      # draw key showing grayscale
      x0 <- 1.5; dx <- 0.15; y0 <- 0; dy <- 0.2; z0 <- 0; dz <- 0.0
      for (i in 0:4) {
        y <- y0 + i*dy
        quads3d( c(x0,x0+dx,x0+dx,x0), c(y,y,y+dy,y+dy), c(z0,z0,z0,z0), col=gray(1-i*0.25) )
      }
      text3d( x0+1.5*dx, y0+dy/2, z0, adj=0, texts="0")
      text3d( x0+1.5*dx, y0+dy*4.5, z0, adj=0, texts=paste( scale ) )
    }
  }
}

invisible( list(counts=tally$counts, nrejects=tally$nrejects, value=value, nties=tally$nties,
           nx=tally$nx, rel.freq=tally$rel.freq, rel.rejects=tally$rel.rejects, mesh=sphere,
           k=k,plot.type=plot.type,positive.only=positive.only,report=report)) }
###############################################################################
histRectangular <- function( x, breaks=10, plot.type="default", freq=TRUE, report="summary", ... ) {
# plot a rectangular grid based histogram for the d-dimensional data in x
# breaks  specifies the bin boundaries in each coordinate, one of:
#  - a vector of length d, which specifies the number of bins in each dimension. Then
#    the bin boundaries in coordinate j are given by seq(min(x[j,]),max(x[j,]),length=breaks[j])
#  - a single number m, which is equivalent to breaks=rep(m,d)
#  - a list with breaks[[j]] a vector of doubles that gives the
#    dividing points/bin boundaries for the j-th coordinate of x

# plot.types to add: flat, sections, none, ..

allowed.types <- c("default", "index", "none", "pillars", "counts" )
stopifnot( is.matrix(x), is.numeric(breaks), is.character(plot.type), is.character(report),
   plot.type %in% allowed.types, report %in% c("all","summary","none") )
n <- nrow(x)
d <- ncol(x)

if( plot.type=="default" ) { plot.type <- ifelse( d <= 3,"pillars","index") }
if( plot.type=="pillars" & d > 3) { stop('plot.type="pillars" does not work for dimension > 3') }
if( plot.type=="counts" & d > 3) { stop('plot.type="counts" does not work for dimension > 3') }
 
# find the range of the data for default bin definition
a <- rep( NA, d )
b <- a
for (j in 1:d) {
  tmp <- range( x[,j] )
  a[j] <- tmp[1]
  b[j] <- tmp[2]
}

# define rectangular bins, convert to H-rep and tally
rect <- RectangularMesh( a, b, breaks=breaks,  )
H <- V2Hrep( rect$S )
tally <- TallyHrep( x, H, report=report )
if (freq) {
  value <- tally$counts
  ylab.string <- "frequency"  
} else {
  value <- tally$relative.freq
  ylab.string <- "relative frequency"   
}

   
# plot results
if ( plot.type=="counts" ) {
  plot( rect, show.labels=TRUE, label.values=value, ... )
} 
if (plot.type=="pillars") {
  if (d==2) {
    open3d()
    DrawPillars( rect$S, value, c(0,0,0), ... ) 
  }
  if (d==3) {
    nbins <- unlist(lapply(rect$breaks,length)) - 1L
    zbreaks <- rect$breaks[[3]]
    z.spacing <-  1.5* max(value)
    z.level <- 0.0
    bins.per.level <- nbins[1]*nbins[2]
    start.S <- 1L
    open3d()
    for (k in 1:nbins[3]) {
      keep <- start.S:(start.S+bins.per.level-1L)
      DrawPillars( rect$S[1:4,1:2,keep], value[keep], c(0,0,z.level), ... ) 
      z.level <- z.level + z.spacing
      start.S <- start.S + bins.per.level
    }
  }
}
if (plot.type=="index") {
  index <- 1:length(value)
  plot(index,value,type='h',ylab=ylab.string,...)
}   

invisible( list(counts=tally$counts, nrejects=tally$nrejects, nties=tally$nties,
          nx=tally$nx, rel.freq=tally$rel.freq, rel.rejects=tally$rel.rejects, mesh=rect,
           plot.type=plot.type,report=report)) }
###############################################################################
histSimplex <- function( x, S, plot.type="default", freq=TRUE, report="summary", ... ) {
# plot a simplex based histogram for the d-dimensional data in x
# S is an array of size (m x d x k), with S[ , ,j] giving the vertices of simplex j

allowed.types <- c("default", "index", "none", "pillars", "counts" )
stopifnot( is.matrix(x), is.numeric(S), is.array(S), is.character(plot.type), is.character(report),
   plot.type %in% allowed.types, report %in% c("all","summary","none")   )
n <- nrow(x)
d <- ncol(x)

if (plot.type == "default") { 
  if (d==2) { plot.type <- "pillars" }
  if (d==3) { plot.type <- "counts" }
  if (d >3) { plot.type <- "index" }
}  
if( plot.type=="pillars" & d != 2) { stop('plot.type="pillars" only works for dimension 2') }
if( plot.type=="counts" & d > 3) { stop('plot.type="counts" does not work for dimension > 3') }
 
# convert from V-representation to H-representation of the simplices and tally
H <- V2Hrep( S )
tally <- TallyHrep( x, H , report=report )
if (freq) {
  value <- tally$counts
  ylab.string <- "frequency"   
} else {
  value <- tally$relative.freq
  ylab.string <- "relative frequency" 
}

# construct "fake" mesh for caller; this is only used internally
vps <- nrow(S)
mvmesh.type <- ifelse( vps==2^d, 7L, 1L ) # hack to detect rectangular mesh
mesh <- list(type="histSimplex",mvmesh.type=mvmesh.type,n=n,vps=vps,S=S)
class(mesh) <- "mvmesh"

# plot results
if (plot.type=="counts") {
  plot( mesh, show.labels=TRUE, label.values=value, ... )
}
if (plot.type=="pillars") {
  open3d()
  DrawPillars( mesh$S, value, ... )
} 
if (plot.type=="index") {
  index <- 1:length(value)
  plot(index,value,type='h',ylab=ylab.string,...)
}   

invisible( list(counts=tally$counts, nrejects=tally$nrejects, nties=tally$nties,
           nx=tally$nx, rel.freq=tally$rel.freq, rel.rejects=tally$rel.rejects, mesh=mesh,
           plot.type=plot.type,report=report)) }
###############################################################################
TallyHrep <- function( x, H, report="summary" ) {
# tally the number of x vectors that are in the simplices specified by H
#
# Input values:
#   x a (d x n) data matrix, columns x[,1],...,x[,n] are data points
#   H a (d x m x nbins) array, H[,,j] is the H-representation of the j-th simplex
#   report should be one of c("summary","all","none") and determines the amount of 
#     reporting for rejects and ties

# input checking
stopifnot( is.matrix(x), is.array(H), length(dim(H))==3, ncol(H)-2==ncol(x),
    report %in% c("all","summary","none")  )

nbins <- dim(H)[3]
bin.counts <- rep( 0L, nbins )
x.counts <- rep( 0L, nrow(x) )
nties <- 0L

# loop through cones, tallying number of data vectors in each cone  
for (i in 1:nbins) {
  satisfy <- SatisfyHrep( x, H[,,i] )
  for (j in satisfy) {
    if (x.counts[j]==0L) {
      bin.counts[i] <- bin.counts[i] + 1L
    } else { 
      # x[j] has already been counted in some cone with index < i
      nties <- nties + 1L          
      if (report=="all") { 
        warning( paste("data vector ",j,"is in multiple cones") )
      }
    } 
    x.counts[j] <- x.counts[j]+1L
  }
}
 
not.in.cone <- which( x.counts==0L ) 
nrejects <- length(not.in.cone) 
if ((nrejects > 0) & (report=="all") ) {
  warning( paste("data values not in any bin: ",paste(not.in.cone,collapse=", ") ) )
}

if ( (report %in% c("summary","all") ) & ((nties > 0) || (nrejects > 0) ) ){
  warning( paste("nties=",nties,"   nrejects=",nrejects) )
}

nx <- nrow(x)
return( list(counts=bin.counts,nrejects=nrejects,nties=nties,nx=nx,
             rel.freq=bin.counts/nx,rel.rejects=nrejects/nx) ) }
###############################################################################
HrepCones <- function( S ) {
# calculate the H-representation for the cones generated by the simplices in S and the origin
# Input:
#   S is a (p x d x ncones) array, with S[,,j] giving the generator of the j-th cone
# Output: 
#   H is a (m x (d+2) x ncones) array, with H[,,k] giving the H-rep. of the k-th cone
#

tmp <- dim( S )
d <- tmp[2]; ncones <- tmp[3]
zeros <- rep(0,d)

for (j in 1:ncones) {
  # set "l" and "b" columns to 0
  hrep <- rcdd::scdd( cbind(  zeros, zeros, S[,,j] ), representation="V" )$output
  if (j==1) {
    H <- array(0.0, dim=c(dim(hrep),ncones) )
  }
  attr( hrep, "representation" ) <- NULL
  H[,,j] <- hrep
}
return( H ) }
###############################################################################
DrawPillars <- function( S, height, shift=rep(0.0,3), ... ) {
# draw "pillars" over the 2d simplices in S
#  S[1:vps,1:2,1:nS] array, with S[,,k] specifying the vertices of the k-th simplex
#  height[1:nS] heights of pillars
#  shift[1:3] shift of the base simplices, typically c(0,0,z0)

if( is.matrix(S) ) { S <- array( S, dim=c(dim(S),1) ) }
stopifnot( is.array(S), length(dim(S))==3, dim(S)[3]==length(height), ncol(S)==2, length(shift)==3 )
vps <- nrow(S); nS <- dim(S)[3]
mvmesh.type <- ifelse( vps==4, 7L, 2L ) # hack so that we know below when there is a rectangle

# loop through simplices, drawing a pillar for each simplex
for (k in 1:nS) {
  bottom <- cbind( S[,,k], rep(0.0,vps)) + matrix( shift, byrow=TRUE, nrow=vps, ncol=3 )
  top <- bottom + matrix( c(0,0,height[k]), byrow=TRUE, nrow=vps, ncol=3 )
  # draw top and bottom
  if (mvmesh.type==7L) {
    # for a rectangular mesh, draw a rectangle; need to reorder vertices
    # warning: this reordering depends on the order in which MeshRectangular defines vertices!
    jj <- c(1,2,4,3,1) 
    lines3d(bottom[jj,1],bottom[jj,2],bottom[jj,3],...)
    lines3d(top[jj,1],top[jj,2],top[jj,3],...)    
  } else { 
    for (i in 1:(vps-1)) {
      for (j in (i+1):vps) {
        lines3d(bottom[c(i, j),1],bottom[c(i, j),2],bottom[c(i, j),3],...) 
        lines3d(top[c(i, j),1],top[c(i, j),2],top[c(i, j),3],...) 
      }
    }
  }   
  
  # connect respective vertices from top to bottom
  for (j in 1:vps) {
    lines3d( c(bottom[j,1],top[j,1]), c(bottom[j,2],top[j,2]), c(bottom[j,3],top[j,3]), ... )
  }
}
aspect3d( 1,1,1 )
}
#########################################################################
histDirectionalQuantileThreshold <- function( x, probs=1, p=2, k=3, positive.only=FALSE, ... ) {
# Plot directional histogram for the data in matrix x with thresholding 
# determined by quantiles of the data.  length(probs) plots are produced, with plot i
# showing the upper extreme probs[i]-percentile of the data.  Thresholding
# is determined using the p-th power norm. k=number of subdivisions.

rowsums <- (rowSums( abs(x)^p ))^(1/p)
thresholds <- quantile( rowsums, 1-probs )
for (i in 1:length(thresholds)) {
  keep <- which( rowsums >= thresholds[i] )
  cat("plot",i,"  threshold=",thresholds[i],"  length(keep)=",length(keep),"   fraction kept=",probs[i],"\n")
  histDirectional( x[keep,], k=k, p=p, positive.only=positive.only, ... )
  title.str <- paste(100*probs[i],"% of data, ",length(keep)," values",sep="" )
  if (ncol(x) == 3) { title3d( title.str ) } else { title( title.str ) }
}}
#########################################################################
histDirectionalAbsoluteThreshold <- function( x, thresholds=0, p=2, k=3, positive.only=FALSE, ... ) {
# Plot directional histogram for the data in matrix x with thresholding 
# determined by abolute threshold.  length(probs) plots are produced, with plot i
# showing the data values .  Thresholding
# is determined using the p-th power norm. k=number of subdivisions.

rowsums <- (rowSums( abs(x)^p ))^(1/p)
for (i in 1:length(thresholds)) {
  keep <- which( rowsums >= thresholds[i] )
  cat("plot",i,"  thresholds=",thresholds[i],"  length(keep)=",length(keep),"\n")  
  histDirectional( x[keep,], k=k, p=p, positive.only=positive.only, ... )
  title.str <- paste("threshold=",thresholds[i],",  ", length(keep)," values",sep="" )
  if (ncol(x) == 3) { title3d( title.str ) } else { title( title.str ) }
}}

