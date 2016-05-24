naSpline<-function(mat,maxgap=1){
  
  ## Make sure input is a matrix
  mat<-as.matrix(mat)
  ## Ensure maxgap is positive
  if(maxgap<1){maxgap<-1;print("maxgap cannot be less than 1. Correcting data using maxgap=1")}
  ## define X range for interpolation
  x<-1:nrow(mat)
  
  ## Find NAs
  na<-is.na(mat)
  ##  Find columns with at least 1 NA
  naCols<-colSums(na)>0
  ## Find how many columns contain NA
  nCols<-sum(naCols)
  ## Exit functin if no NAs are found
  if(nCols==0){print("No data needs correcting");return(mat)}
  
  ## exclude columns that have no NAs
  smallNa<-na[,naCols,drop=FALSE]
  ## exclude columns that have no NAs
  smallMat<-mat[,naCols,drop=FALSE]
  
  ## Loop over these columns and...
  for(i in 1:nCols){
    ## Identify runs of NAs
    rl<-rle(smallNa[,i])
    ## only select runs less than maxgap
    rl$values<-rl$values & rl$lengths <=maxgap
    ## Identify indices of these runs
    ind<-inverse.rle(rl)
    ## Interpolate at these indices... if there is something to interpolate
    if(sum(ind)>0){
    smallMat[ind,i]<-spline(x = x,y = smallMat[,i],xout = x[ind],n = length(x))$y
    }
    }
  ## Put the interpolated subset back into the original matrix
  mat[,naCols]<-smallMat
  ## return corrected matrix
  return(mat)
}
