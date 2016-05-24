#' @title Non cumulative similarity function
#' @description 
#' To calculate non cumulative similarity among common and group loadings
#' @param x a numeric vector
#' @export
#' @keywords internal
similarity_noncum<-function(loadings_matrices, NAMES){
  
  nb= length(loadings_matrices)
  H=ncol(loadings_matrices[[1]])
  sim=vector("list",H)
  
  for(h in 1:H){
    MM=matrix(1, nrow=nb, ncol=nb)
    for(aa in 1:(nb-1)){
      for(bb in (aa+1):nb){
        cc= abs(as.numeric(t(loadings_matrices[[aa]][,h]) %*% loadings_matrices[[bb]][,h]))
        MM[aa,bb]=cc
        MM[bb,aa]=MM[aa,bb]
      }
    }
    MM=MM
    
    
    colnames(MM)= NAMES
    rownames(MM)= NAMES
    
    
    sim[[h]]=round(MM,3)
  }
  return(sim)
}
