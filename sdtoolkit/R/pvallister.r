
#nullprob wrapper - applies null prob to each of several boxes in a trajectory


pvallister <- function(checpts,x,y,trajinf){

  nboxes <- length(checpts)
  pvallist <- list()
  ofboxlist <- list()
  
  #first convert all boxes to old form:
  
  for (boxind in checpts){
    
    boxy <- list(box=trajinf$box[[boxind]],dimlist=trajinf$dimlist[[boxind]])
    ofboxlist[[boxind]] <- boxconverter(boxy)
    
  }
  
  for (i in checpts){
  
    pvallist[[i]] <- cbind(ofboxlist[[i]][[1]],nullprob(cbind(x,y),y=NULL,ofboxlist[[i]]))
  
  }
  
  return(pvallist)

}

