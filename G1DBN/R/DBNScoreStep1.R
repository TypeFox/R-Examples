## __________________________________________________________
##
## FUNCTION DBNScoreStep1
##
## Given a time series dataset for p genes, this function infers a 1st-order 
## dependence score matrix S1 (p x p) which contains the score of each edge 
## of a Dynamic Bayesian Network (DAG G(1)) describing first order dependencies 
## between successives variables. 
##
## The smallest score points out the most significant edge for the 1st order 
## dependence DAG G(1). 
## 
## The sets of both predictor and target genes can be reduced to different 
## subsets of the p genes. 
##
## DBNScoreStep1 is the first step of the estimation procedure described in the
## references. See function DBNScoreStep2 to perform the second step selection 
## and infer a score matrix describing full order dependencies.
## __________________________________________________________
##


DBNScoreStep1<-function(data, method='ls',predPosition=NULL,
                        targetPosition=NULL) {

  ## ===============================================
  ## INITIALIZING
  ## _______________________________________________

  data<-as.matrix(data)
  n=dim(data)[1] # nb of time points
  p=dim(data)[2] # nb of genes
  
  ## If predictor or target genes position is specified,
  ## we just work of them
  if(is.null(predPosition)){
    pred <- 1:p}
  else {
    pred <- predPosition}
  d <- length(pred)
  
  if(is.null(targetPosition)){
    tar <- 1:p}
  else {
    tar <- targetPosition}
  r <- length(tar)

  ## The score matrix
  S1ls=NULL
  S1tukey=NULL
  S1huber=NULL
  
  if('ls' %in% method){
    S1ls<-matrix(0,r,d)
  }
  if('tukey' %in% method){
    S1tukey<-matrix(0,r,d)
  }
  if('huber' %in% method){
    S1huber<-matrix(0,r,d)
  }

  ## ===============================================
  ## BUILDING THE SCORE MATRICES
  ## _______________________________________________

  ## Print the total number of vertices

  cat("Treating", r ,"vertices:\n")
  cpt=10
 
  for (i in 1:r){
    ## Print percentage of accomplished vertices
     if ( ((i/r)*100)>=cpt ) {
      cat(cpt, "% ",sep="")
      cpt=cpt+10
    }
  
   
    ## The regression model is
    ## Y = X
    
    ## Y is the vector containing the target genes
    ## on which the regression will be performed
    ## The time point 1 is removed
    y<-data[2:n,tar[i]]

    for (j in  1:(d-1)){
      
      ## for all k > j
      for (k in c(1:d)[-c(1:j)]){

        ## X is a matrix with two columns containing the
        ## predicted gene on which the regression will
        ## be performed. The two columns contain
        ## respectively the data corresponding to the jth
        ## and the kth tested gene.
        ## The time n is removed
        x<-data[1:(n-1),c(pred[j],pred[k])]
        
        ## =====================================
        ## ESTIMATION...
        ## 
        ## Three estimators are available :
        ## - Least square
        ## - Huber
        ## - Tuckey
        
        ## =====================================
        ## LEAST SQUARE ESTIMATOR
        if('ls' %in% method){
          lm.3<-lm(y~x)
          prob<-abs(summary(lm.3)$coef[,"Pr(>|t|)"])  
          ## coefficient aij(k) : aij given k 
          S1ls[i,j]<-max(prob[2],S1ls[i,j],na.rm=TRUE)
          ## coefficient aik(j) : aik given j 
          S1ls[i,k]<-max(prob[3],S1ls[i,k],na.rm=TRUE)
        }
        
        ## =====================================
        ## TUKEY'S ESTIMATOR
        if('tukey' %in% method){
          bisq.3<-rlm(y~x,method='MM')
          prob<-2*(1-pt(abs(summary(bisq.3)$coef[,"t value"]),n-4))
          ## coefficient aij(k) : aij given k 
          S1tukey[i,j]<-max(prob[2],S1tukey[i,j],na.rm=TRUE)
          ## coefficient aik(j) : aik given j 
          S1tukey[i,k]<-max(prob[3],S1tukey[i,k],na.rm=TRUE)
        }

        ## =====================================
        ## HUBER'S ESTIMATOR
        if('huber' %in% method ){
          hub.3<-rlm(y~x)
          prob<-2*(1-pt(abs(summary(hub.3)$coef[,"t value"]),n-4))
          ## coefficient aij(k) : aij given k 
          S1huber[i,j]<-max(prob[2],S1huber[i,j],na.rm=TRUE)
          ## coefficient aik(j) : aik given j 
          S1huber[i,k]<-max(prob[3],S1huber[i,k],na.rm=TRUE)
       }

         
      } # end k
    } # end j
  } # end i 
  cat("\n")

  ## The score matrices are return
  list(S1ls=S1ls,S1tukey=S1tukey,S1huber=S1huber)
}
