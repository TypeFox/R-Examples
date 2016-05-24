.ulapply <- function( l, f, ... ) { 
  unlist(lapply(l,f, ...), use.names=FALSE) 
}

.setColAlpha <- function( col, alpha=128 ){
	do.call( rgb, c(split(col2rgb(col),c("red","green","blue")),
		alpha=alpha,maxColorValue=255) )
}

.minMaxOfTwoMatrices <- function(a, b, func) {
  #   print("a:")
  #   print(a)
  #   print("b:")
  #   print(b)
  res <- matrix(nrow=nrow(a), ncol=ncol(a))
  #   print(ncol(a))
  for (j in 1:ncol(a)) {
    #     print(j)
    res[1, j] <- min(a[1, j], b[1, j])
    res[2, j] <- max(a[2, j], b[2, j])
  }
  colnames(res) <- colnames(a)
  rownames(res) <- rownames(a)
  return(res)
}

.boundingBoxTrack <- function(track) {
  rbind(min = apply(track, 2, min), max = apply(track, 2, max))
}

.subtrackIndices <- function( x, overlap, lengths=nrow(x) ){
  if (overlap > x-1) {
    warning("Overlap exceeds segment length")
    overlap <- x-1
  }
  xl <- c(0,cumsum(lengths))
  r <- .ulapply(seq_along(lengths),function(l){
  	if( xl[l+1]-xl[l] <= x ){
  		c()
  	} else {
  		seq(xl[l]+1,xl[l+1]-x,x-overlap)
  	}
  })
  cbind( first=r, last=r+x )
}


.computeSegmentwiseMeans <- function(track, measure, min.segments=1, only.for.i=NULL, ...) { 
  if (nrow(track) <= min.segments) {
    return(NA)
  }
  means <- c()
  debug.nr.tracks <- 0
  if (is.null(only.for.i)) {
    for (i in (min.segments):(nrow(track) - 1)) {
      subtracks <- .computeSubtracksOfISegments(track, i, ...)
      val <- sapply(subtracks, measure)
      means[i - min.segments + 1] <- mean(val)#, na.omit=TRUE)    # na.rm produces NaN if all values are NA
    }
  } else {
    if ((only.for.i >= min.segments) && (only.for.i <= (nrow(track) - 1))) {
      subtracks <- .computeSubtracksOfISegments(track, only.for.i)
      val <- sapply(subtracks, measure)
      means <- mean(val)#, na.omit=TRUE)
    } else {
      return(NA)
    }
  }
  return(means)
}

.computeSubtracksOfISegments <- function(track, i, overlap=i-1) {
  if (nrow(track) <= i) {
    return( as.tracks(list()) )
  } 
  if (overlap > i-1) {
    warning("Overlap exceeds segment length")
    overlap <- i-1
  }
  l <- list()
  for (j in seq(1, (nrow(track)-i), max(1,i-overlap))) {
    l[[as.character(j)]] <- track[j:(j+i), ,drop=FALSE]
  }    
  return( as.tracks(l) )
}


.segmentwiseMeans <- function(measure, ...) {
  return(function(track) {
    .computeSegmentwiseMeans(track, measure, ...)
  })
}

# Helper function of beauWalker() and is not directly called by user
# This function samples uniformly from unit sphere
.beaucheminSphereSampler <- function(d=3){
  x <- rnorm(d)
  return(x/sqrt(sum(x^2)))
}

# Helper function of beauWalker() and is not directly called by user
# This function samples from the trianglular distribution
.beaucheminTriangularSampler <- function() sqrt(4*runif(1))-1

# Creates a matrix that will rotate unit vector a onto unit vector b
.beaucheminRotationMatrix <- function(a,b){
		# handle undefined cases
	    theta <- acos( sum(a*b) )
    	R <- diag(rep(1,3))
		if( theta < 0.001 ){
			return(R)
		}
		if( pi-theta < 0.001 ){
			return(-R)
		}
		# compute normalized cross product of a and b
		x <- c(a[2]*b[3]-a[3]*b[2],
			a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
		x <- x / sqrt(sum(x^2))
		A <- rbind( c(0,-x[3],x[2]),c(x[3],0,-x[1]),c(-x[2],x[1],0) )
		R <- R + sin(theta)*A + (1-cos(theta))* (A%*%A)
    	return(R)
}

# Helper function of beauWalker() and is not directly called by user
# returns a direction for the cell to travel based on model parameters specified in beauWalker()
.beaucheminPickNewDirection <- function(
	old.direction,
	p.bias,p.persist,bias.dir,taxis.mode,
	t.free,v.free,rot.mat){
	if( runif(1) < p.persist ){
		return(c(old.direction, t.free))
	}
	d <- .beaucheminSphereSampler(3)
	if( taxis.mode == 0 ){
		return(c(d,t.free))
	} else if(taxis.mode==1){
		# orthotaxis
		return(c(v.free * d*(1+p.bias*sum(d*bias.dir)),t.free))
	} else if(taxis.mode == 2){
		# topotaxis	
		if( runif(1) < p.bias ){
			# Approach: generate new direction as if the bias direction were (1,0,0).
			# Then rotate the resulting direction by the angles between (1,0,0) and the true
			# bias direction. 
			d[1] <- .beaucheminTriangularSampler()
			circ <- .beaucheminSphereSampler(2)
			circle.scale <- sqrt(1-d[1]^2)
			d[2:3] <- circ[1:2]*circle.scale
			return(c( v.free * rot.mat%*%d, t.free))
		} else {
			return(c( v.free * d, t.free))
		}
	} else if(taxis.mode == 3){
	 	 # klinotaxis
		return(c(d*v.free,t.free*(1+p.bias*sum(d*bias.dir))))
	}
}

.gaps <- function(x, tol=0.05, deltaT=NULL){
	if( nrow(x) < 2 ){
		return(integer())
	}
	gapThreshold <- tol*deltaT
	xt <- x[,1]
	xt[apply( is.na(x), 1, any )] <-  NA
	r <- abs(diff(x[,1])-deltaT) > gapThreshold
	r[is.na(r)] <- TRUE
	which(r)
}

.normalizeTrack <- function(track) {
	cbind( track[,1,drop=FALSE], sweep(track[,-1],2,track[1,-1]) )
}

