#' Title
#'
#' @param nrows 
#' @param ncols 
#' @param rook 
#'
#' @return
#' @export
#'
#' @examples
lat2w <- function(nrows=5, ncols=5, rook=TRUE){
  
  n <- nrows*ncols
  
  m <- matrix(1:n,nrows,ncols, byrow=TRUE)
  
  w <- matrix(0,n,n)
  
  x <- vector(mode="numeric", length=0)
  
  nbs <- rep(list(x), n)
  
  counter<-1
  
  for(i in 1:nrows){
    for(j in 1:ncols){
      nb <- vector(mode="numeric", length=0)
      if((j-1)>0) {
        nb<-c(nb,m[i,j-1])
        w[counter,m[i,j-1]]<-1
        } #East
      if((j+1)<=ncols) {
        nb<-c(nb,m[i,j+1])
        w[counter,m[i,j+1]]<-1
        } #West
      if((i-1)>0) {
        nb<-c(nb,m[i-1,j])
        w[counter,m[i-1,j]]<-1
        } #North
      if((i+1)<=nrows) {
        nb<-c(nb,m[i+1,j])
        w[counter,m[i+1,j]]<-1
        } #South
      if(rook==FALSE)
        {
        if((i-1)>0 && (j-1)>0) {
          nb<-c(nb,m[i-1,j-1])
          w[counter,m[i-1,j-1]]<-1
        } #NorthEast
        if((i-1)>0 && (j+1)<=ncols) {
          nb<-c(nb,m[i-1,j+1])
          w[counter,m[i-1,j+1]]<-1
        } #NorthWest
        if((i+1)<=nrows && (j-1)>0) {
          nb<-c(nb,m[i+1,j-1])
          w[counter,m[i+1,j-1]]<-1
        } #SowthEast
        if(i+1<=nrows && (j+1)<=ncols) {
          nb<-c(nb,m[i+1,j+1])
          w[counter,m[i+1,j+1]]<-1
        } #SouthWest
      }
    nbs[[counter]]<-nb
    counter<-counter+1
    }
  }
  return(list(nbs=nbs,w=w))
}