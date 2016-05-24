## __________________________________________________________
##
## FUNCTION BuildEdges
##
## Given a score matrix, this function builds the list of the
## edges of the associated network. Edges are ordered according 
## to their scores.
##
## The score matrix has been computed from a network inference 
## algorithm (e.g. DBNScoreStep1 or DBNScoreStep2, Shrinkage, 
## Lasso, ...). 
##
## An optional threshold can be specified, as well as a maximal
## number of edges.
##
##
## __________________________________________________________
##
BuildEdges <- function(score,threshold=1,nb=NULL,
  targetNames=NULL,predNames=NULL,prec=3,dec=FALSE){

  ## ===============================================
  ## INITIALIZING
  ## _______________________________________________
  r <- dim(score)[1] # nb of target genes 
  d <- dim(score)[2] # nb of predictor genes
  
  ## ===============================================
  ## BUILDING THE MATRIX OF EDGES
  ## _______________________________________________  
  
  ## ordering the scores
  scoreIdx <- order(score,decreasing=dec)
  nbEdges <- min(nb,sum((score<threshold)*(dec==FALSE)+(score>threshold)*(dec==TRUE),na.rm=TRUE),length(scoreIdx))
  
  ## coordinate in the matrix
  Pred <- ceiling(scoreIdx[1:nbEdges]/r) 
  Target <- scoreIdx[1:nbEdges]-(Pred-1)*r
 
  ## names if specified
  if(!(is.null(predNames))){ 
    Pred <- predNames[Pred]
  }
  if(!(is.null(targetNames))){ 
    Target <- targetNames[Target]
  }
   
  ## Returning the matrix of Edges
  Score=round(score[scoreIdx][1:nbEdges],prec)
  return(as.matrix(cbind(Pred,Target,Score)))
}
