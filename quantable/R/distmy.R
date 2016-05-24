#' distance among matrix columns with freely choosable distance function
#' @param x data
#' @param func function taking 2 arrays x, y
#' @param init how to initialize the output matrix
#' @param diag should the diagonal be also computed
#' @return matrix with dist[i,j] = func(x[,i],[x,j])
#' @export
#' @examples
#' mat = matrix(rnorm(10*5000),ncol=10)
#' redist = distmy(mat,function(x,y){mean(abs(x-y))},init=NA,diag=FALSE)
#' image(redist)
#' redist = distmy(mat,cor,init=0,diag=FALSE)
#' image(redist)
#' redist = distmy(mat,function(x,y){ks.test(x,y)$p.value},init=1,diag=TRUE)
#' image(redist)
#' hist(uppertriang(redist))
#' range(uppertriang(redist))
#' which(redist < 0.05 , arr.ind = TRUE)
distmy = function( x, func, init=NA , diag = TRUE)
{
  f = NULL
  if(diag){
    f = function(i,j){i<=j}
  }else{
    f = function(i,j){i<j}
  }
  nout = dim(x)[2]
  resP = matrix(init,nrow=nout,ncol=nout)
  for(i in 1:nout){
    for(j in 1:nout){
      if(f(i,j)){
        test = func(x[,i],x[,j])
        resP[i,j] = test
      }
    }
  }
  resP[lower.tri(resP)] = t(resP)[lower.tri(resP)]
  rownames(resP) = colnames(x)
  colnames(resP) = colnames(x)
  return(resP)
}

