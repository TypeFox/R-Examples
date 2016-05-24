




LKrig.normalize.basis.slow <- function(Level, LKinfo, PHItemp){
  # get  SAR matrix at level Level
   tempB <- LKrig.MRF.precision(LKinfo$mx[Level], LKinfo$my[Level], 
                     a.wght = (LKinfo$a.wght)[[Level]], edge = LKinfo$edge, 
                     distance.type = LKinfo$distance.type)   
               tempB <- LKrig.spind2spam(tempB)
# quadratic form involves applying inverse precision matrix to basis function evaluted at
# locations for evaluation
               wght <- LKrig.quadraticform( t(tempB) %*% tempB, PHItemp)
return( wght)}


LKrig.quadraticform <- function(Q, PHI) {
    #   finds     the quadratic forms    PHI_j^T Q.inverse PHI_j  where PHI_j is the jth row of
    #   PHI.
    #   The loop is to avoid using memory for the entire problem if more than 2000 elements.
    nrow <- dim(PHI)[1]
    choleskyPrecision<- chol(Q)
    if (nrow > 1) {
        BLOCKSIZE <- 2000
        wght <- rep(NA, nrow)
        counter <- 1
        while (counter < nrow) {
            ind <- counter:min((counter + BLOCKSIZE), nrow)
            A <- forwardsolve(l = choleskyPrecision, transpose = TRUE, 
                x = t(PHI[ind, ]), upper.tri = TRUE)
            wght[ind] <- c(colSums(A^2))
            counter <- counter + BLOCKSIZE
        }
    }
    else {
        # Unfortunately the case when there is one row needs to be handled separately.
        A <- forwardsolve(l = choleskyPrecision, transpose = TRUE, x = t(PHI), upper.tri = TRUE)
        wght <- sum(A^2)
    }
    return(wght)
}




LKrig.normalize.basis.fast<- function(Level, LKinfo, x ){
  if( !(LKinfo$fastNormalization)){
     stop("Not setup for fast normalization")}
  # convert the locations to the integer scale of the lattice
  # at the level == Level
  mxLevel<- LKinfo$mx[Level]
  myLevel<- LKinfo$my[Level]
  
  gridStuff<- LKinfo$grid[[Level]]
  xmin<- gridStuff$x[1]
  ymin<- gridStuff$y[1]
  dx<-  gridStuff$x[2]-  gridStuff$x[1]
  dy<-  gridStuff$y[2]-  gridStuff$y[1]            
  xLocation<- scale( x, center= c( xmin, ymin), scale= c( dx, dy)) + 1
  nLocation<- nrow( xLocation)
  setupList<- LKinfo$NormalizeList[[Level]]
  return(
         .Fortran("findNorm",
                          mx = as.integer(mxLevel),
			  my = as.integer(myLevel),
		      offset = as.double(LKinfo$overlap),
			  Ux = as.double(setupList$Ux),
			  Dx = as.double(setupList$Dx),
			  Uy = as.double(setupList$Uy),
			  Dy = as.double(setupList$Dy),
                   nLocation = as.integer(nLocation),
                   xLocation = as.double( xLocation),
                     weights = as.double( rep(-1,nLocation) ), 
	       	           Z = matrix(as.double(0),mxLevel,myLevel)
                  )$weights
         )
}

LKrig.make.Normalization<- function(mx,my, a.wght){
  out<- list()
  nlevel<- length( a.wght)
  for ( l in 1:nlevel){
    out<- c(out, list(LKrigMRFDecomposition( mx[l], my[l], a.wght[[l]] )) )
  }
  return(out)
 } 

LKrigMRFDecomposition<-  function( mx,my,a.wght){
  Ax<- diag( a.wght/2, mx)
  Ax[  cbind( 2:mx, 1:(mx-1)) ] <- -1
  Ax[  cbind( 1:(mx-1), 2:mx) ] <- -1
# 
  Ay<- diag( a.wght/2, my)
  Ay[  cbind( 2:my, 1:(my-1)) ] <- -1
  Ay[  cbind( 1:(my-1), 2:my) ] <- -1
  eigen( Ax, symmetric=TRUE) -> hold
  Ux<- hold$vectors
  Dx<- hold$values
  eigen( Ay, symmetric=TRUE) -> hold
  Uy<- hold$vectors
  Dy<- hold$values
#
  return( list( Ux=Ux, Uy=Uy, Dx=Dx, Dy=Dy ))
       }

LKrig.normalize.basis <- function(Level, LKinfo, PHI){
            tempB <- LKrig.MRF.precision(LKinfo$mx[Level], LKinfo$my[Level], 
                               a.wght = (LKinfo$a.wght)[[Level]],
                               edge = LKinfo$edge, 
                               distance.type = LKinfo$distance.type)
            tempB <- LKrig.spind2spam(tempB)
            Q.level <- t(tempB) %*% tempB
            wght <- LKrig.quadraticform(Q.level, PHI)
            return(wght)
}
