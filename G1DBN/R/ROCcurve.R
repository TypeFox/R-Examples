## __________________________________________________________
##
## FUNCTION ROCcurve
##
## Given a score matrix and a validation matrix, this function
## allows to compute the corresponding ROC curve by returning a list
## with the respective x and y coordinates of the ROC curve. 
##
## The score matrix has been computed from a network inference algorithm 
## (e.g. DBNScoreStep1 or DBNScoreStep2, Shrinkage, Lasso, ...).
## __________________________________________________________
##

ROCcurve <- function(score,validMat,dec=FALSE){

 ## Initializing...
 r <- dim(score)[1] # nb of target genes
 d <- dim(score)[2] # nb of predictor genes

 ## TestEdges is a boolean vector that contains
 ## the comparison between inferred and real edges
 ScoreIdx <- order(score,decreasing=dec)

 TrueEdges <- which((abs(validMat)>0))
 TestEdges <- ScoreIdx %in% TrueEdges

 ## Building the ROC curve...
 rocx<-matrix(0,length(ScoreIdx)+1,1) # x coord of the roc curve
 rocy<-matrix(0,length(ScoreIdx)+1,1) # y coord of the roc curve

 for (i in 1:length(ScoreIdx)){
   if(TestEdges[i]==TRUE){
     rocx[i+1]<-rocx[i]
     rocy[i+1]<-rocy[i]+1}
   else{
     rocx[i+1]<-rocx[i]+1
     rocy[i+1]<-rocy[i]
   }
 }

 return(list(x=rocx,y=rocy))
}

