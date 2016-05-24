## __________________________________________________________
##
## FUNCTION PRCurve
##
## Given a score matrix and a validation matrix, this function
## allows to compute the corresponding Precision-Recall (PR) curve by
## returning a list with the respective x coordinates (recall) and
## y coordinates (precision) of the PR curve. The recall is equal to
## the sensitivity, that is the number of true positive out of the number
## of true edges to de detected. The precision is the Positive Predictive
## Value, that is the number of true positive edges out of the number of 
## selected edges. 
##
## The score matrix has been computed from a network inference algorithm 
## (e.g. DBNScoreStep1 or DBNScoreStep2, Shrinkage, Lasso, ...).
## __________________________________________________________
##

PRcurve<-function(score,validMat,dec=FALSE){

  ## Initializing...
  p<-dim(score)[1]
  q<-dim(score)[2]
  nbPos=sum(validMat)
  
  ## Odering the scores
  tri<-sort(score,decreasing=dec)

  precision=array(0,length(tri))
  recall=array(0,length(tri))
  
  ## Building the PR curve...
  for(j in 1:length(tri)){
    if(length(tri)>10){
      if ((j %% (round(length(tri)/10))) == 0) {
        cat(round(10*j/(round(length(tri)/10))),"% ")
      }
    }
    recall[j]=sum(validMat*(score<=tri[j]), na.rm = TRUE)/nbPos
    precision[j]=sum(validMat*(score<=tri[j]), na.rm = TRUE)/sum(score<=tri[j] *(score<1), na.rm = TRUE)
  }
  cat("\n")
  return(list(recall=recall,precision=precision))
}

