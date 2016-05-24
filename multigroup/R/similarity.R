#' @title similarity function
#' @description 
#' To calculate similarity among common and group loadings
#' @param x a numeric vector
#' @export
#' @keywords internal
similarity_function<-function(loadings_matrices, NAMES){
  
  nb= length(loadings_matrices)
  H=ncol(loadings_matrices[[1]])
  sim=vector("list",H)
  
  for(h in 1:H){
    MM=matrix(h, nrow=nb, ncol=nb)
    for(aa in 1:(nb-1)){
      for(bb in (aa+1):nb){
        cc=0
        for(i in 1:h){
          cc=cc+abs(as.numeric(t(loadings_matrices[[aa]][,i]) %*% loadings_matrices[[bb]][,i]))
        }
        MM[aa,bb] = cc
        MM[bb,aa] = MM[aa,bb]
      }
    }
    MM = MM/h
    
    
    colnames(MM)= NAMES
    rownames(MM)= NAMES
    
    
    sim[[h]]=round(MM,3)
  }
  return(sim)
}
