## __________________________________________________________
##
## FUNCTION DBNScoreStep2
##
## Given a time series dataset for p genes, a 1st order dependence 
## score matrix S1 (obtained with function DBNScoreStep1) and a 
## threshold alpha1 for edge selection in matrix S1, this function 
## infers the score of each edge of a Dynamic Bayesian Network (DAG G) 
## describing full order dependencies between successive variables. 
## 
## This is the second step of the inference procedure described in 
## the references. 1st step DBNScoreStep1 allows to reduce the number 
## of potential edges, DBNScoreStep2 performs the last step selection. 
## 
## The smallest score points out the most significant edge. 
## __________________________________________________________
##


DBNScoreStep2<-function(S1,data,method='ls',alpha1,
                        predPosition=NULL,targetPosition=NULL){

  ## ===============================================
  ## INITIALIZING
  ## _______________________________________________

  S1<-as.matrix(S1)
  data<-as.matrix(data)
  r <- dim(S1)[1] # nb of target genes 
  d <- dim(S1)[2] # nb of predictor genes
  
  if((r!=dim(data)[2]) && is.null(targetPosition)) {
    stop("Please specify the target genes position in the dataset")
  }
  
  if((d!=dim(data)[2])&& is.null(predPosition)){
    stop("Please specify the predictor genes position in the dataset")
  }
  
  ## If predPosition and / or target genes are specified,
  ## we just work of them

  ## vector of predictor genes
  if(is.null(predPosition)){
    pred <- 1:d}
  else {
    pred <- predPosition}
  
  ## vector of target genes 
  if(is.null(targetPosition)){
    tar <- 1:r}
  else {
    tar <- targetPosition}
  
  n <- dim(data)[1] # nb of time points
  S2 <- matrix(NA,r,d) # The score matrix - 2nd step
  
  ## ===============================================
  ## BUILDING THE SCORE MATRIX
  ## _______________________________________________
  
  ## nb of parents for each target gene in G1
  NbPar <- apply(S1<alpha1,1,sum) 
  
  if (max(NbPar)>(n-1)){
    print(paste("Warning: Threshold alpha1 is two high, some nodes have more than ",n-1," parents in the inferred DAG G(1)"))
  } else {
    
    ## for genes having parents in G1
    for (i in which(NbPar>=1)){
      
      ## Considering only predictors such as S1 < alpha1
      selec1<-which(S1[i,]<alpha1)
      
      ## The regression model is
      ## Y = X
           
      ## Y is the vector containing the target genes
      ## on which the regression will be performed
      ## The time point 1 is removed
      y <- data[2:n,tar[i]]
      
      ## X is a matrix with two columns containing the
      ## predicted gene on which the regression will
      ## be performed.
      ## The time n is removed
      x <- matrix(data[1:(n-1),pred[selec1]],n-1,length(selec1))

      ## The three following estimators are available :
      ## - Least square
      ## - Huber
      ## - Tuckey
      
      ## =====================================
      ## LEAST SQUARE ESTIMATOR
      if ('ls' %in% method) {
        ls<-lm(y~x)
        if(length(summary(ls)$coeff[,"Pr(>|t|)"])!=(length(selec1)+1)){
          print(paste("pb regression gene",i))
        } else {
          S2[i,selec1]<-summary(ls)$coeff[2:(length(selec1)+1),"Pr(>|t|)"]
        }
      }

      ## =====================================
      ## TUKEY'S ESTIMATOR
      if ('tukey' %in% method) {
        tuk<-rlm(y~x,method='MM')
        if(length(summary(tuk)$coef[,"t value"])!=(length(selec1)+1)){
          print(paste("pb regression gene",i))
        } else {
          S2[i,selec1]<-2*(1-pt(abs(summary(tuk)$coef[2:(length(selec1)+1),"t value"]),n-length(selec1)-1))
        }
      }
      
      ## =====================================
      ## HUBER'S ESTIMATOR
      if ('huber' %in% method) {
        hub<-rlm(y~x)
        if(length(summary(hub)$coef[,"t value"])!=(length(selec1)+1)){
          print(paste("pb regression gene",i))
        } else {
          S2[i,selec1]<-2*(1-pt(abs(summary(hub)$coef[2:(length(selec1)+1),"t value"]),n-length(selec1)-1))
        }
      }
      
    } # end for i
  } # end if/else
  
  ## The score matrix is return
  return(S2)
}
