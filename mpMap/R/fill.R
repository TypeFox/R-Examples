fill <-
function(mat, missfx=2, ...)
{
  mat2 <- mat
  upperTriangle(mat2) <- -1

  ## basically want to find the closest entry by index in the matrix which
  ## is not missing and fill in - note that if the

  if (missfx==1) {
    missa <- which(is.na(mat2)) 
    indexa <- vector(length=length(missa))
    for (i in 1:length(missa))
    indexa[i] <- which(!is.na(mat2))[which.min(abs(which(!is.na(mat2))-missa[i]))]
    mat2[missa] <- mat2[indexa]
  }

  ### fill markers in with the average of the distance to two markers away 
  ### - the distance between the two further away markers. (note, this assumes
  ### that there are no other missing values nearby...)
  if (missfx==2) {
    mat2[mat2==.5] <- .49
    missa <- which(is.na(mat2), arr.ind=T)
    if (nrow(missa)>0)
    for (i in 1:nrow(missa)) {
	dist1 <- dist2 <- NA

	### can we take care of the situation where those next points are missing too? 
	k=0
	while((missa[i,1]+k)<ncol(mat2) & is.na(dist1)) {
	dist1 <- haldaneX2R((haldaneR2X(mat2[missa[i,1]+k, missa[i,2]])-haldaneR2X(mat2[missa[i,1]+k, missa[i,1]]))) 
	k <- k+1
	}
	if (!is.na(dist1) & dist1<0) dist1 <- NA
	p=0
	while((missa[i,2]-p)>1 & is.na(dist2)) {
	dist2 <- haldaneX2R((haldaneR2X(mat2[missa[i,1], missa[i,2]-p])-haldaneR2X(mat2[missa[i,2], missa[i,2]-p])))
	p <- p+1
	}
	if (!is.na(dist2) & dist2<0) dist2 <- NA
	if (is.na(dist1)) mat2[missa[i,1], missa[i,2]] <- dist2
	if (is.na(dist2)) mat2[missa[i,1], missa[i,2]] <- dist1
	if (!is.na(dist1) & !is.na(dist2)) {
		if (p>k) mat2[missa[i,1], missa[i,2]] <- dist1
		if (p<k) mat2[missa[i,1], missa[i,2]] <- dist2
		if (p==k) mat2[missa[i,1], missa[i,2]] <- (dist1+dist2)/2
	}
    }
  }

  lowerTriangle(mat) <- lowerTriangle(mat2)
  upperTriangle(mat) <- upperTriangle(t(mat2))
  
  return(mat)
}
