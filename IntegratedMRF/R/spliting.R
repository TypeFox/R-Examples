#' Split of the Parent node
#' 
#' Split of the training samples of the parent node into the child nodes based on the feature and threshold that produces the minimum cost
#'  
#' @param X Input Training matrix of size M x N, M is the number of training samples and N is the number of features
#' @param Y Output Training response of size M x T, M is the number of samples and T is the number of output responses
#' @param m_feature Number of randomly selected features considered for a split in each regression tree node.
#' @param Index Index of training samples
#' @param Cov_Y Covariance matrix of Output Response matrix for MRF(Give Zero for RF)
#' @param Command 1 for univariate Regression Tree (corresponding to RF) and 2 for Multivariate Regression Tree (corresponding to MRF)
#' @return List with the following components:
#' \item{index_left}{Index of the samples that are in the left node after splitting}
#' \item{index_right}{Index of the samples that are in the right node after splitting}
#' \item{which_feature}{The number of the feature that produces the minimum splitting cost}
#' \item{threshold_feature}{The threshold value for the node split. 
#' A feature value less than or equal to the threshold will go to the left node and it will go to the right node otherwise.}
#' @details 
#' At each node of a regression a tree, a fixed number of features (m_feature) are selected randomly to be 
#' considered for generating the split. Node cost for all selected features along with possible n-1 thresholds for 
#' n samples are considered to pick the feature and threshold with minimum cost.   
#' @examples
#' X=matrix(runif(20*100),20,100)
#' Y=matrix(runif(20*3),20,3)
#' m_feature=5
#' Index=1:20
#' Cov_Y=stats::cov(Y)
#' Command=2#MRF, as number of output feature is greater than 1
#' Split_criteria=spliting(X,Y,m_feature,Index,Cov_Y,Command) 
#' @export
spliting <- function(X,Y,m_feature,Index,Cov_Y,Command){
  x=X[Index, ]
  if (Command==1){
    y=matrix(Y[Index, ],ncol=1)
  }else {
    y=Y[Index, ]
  }
  f = ncol(x) # number of features
  f_1 =sort(sample(f, m_feature)) #randomly taken 10 features, for each splits vary
  N = nrow(x)
  min_score =NULL
  DL=matrix(data=NA, nrow=(N-1), ncol=length(f_1))
  DR=matrix(data=NA, nrow=(N-1), ncol=length(f_1))
  id =NULL
  for(j in 1:length(f_1)){
    xj = x[ ,f_1[j]]
    tmp=sort(xj,index.return=TRUE)   # sort the values
    s=tmp$x
    idx=tmp$ix 
    id2=NULL
    for (k1 in 1:(N-1)){
      k=k1
      yk_left=matrix(y[idx[1:k1],],ncol=ncol(y)) # will create problem with mrf
      yk_right=matrix(y[idx[(k1+1):length(idx)],],ncol=ncol(y))
      D_left= Multi_D_mod(yk_left,Cov_Y,Command)
      D_right= Multi_D_mod(yk_right,Cov_Y,Command)
      D=D_left+D_right
      DL[k,j]=D_left
      DR[k,j]=D_right
      if(j==1 && k==1){
        min_score = D
        which_feature = f_1[1]
        indexX=c(k,j)
        threshold_feature = (s[k] + s[k+1])/2
        ij_i = 1
        ij_j = 1
      }
      
      if(D< min_score){
        min_score = D
        which_feature = f_1[j]
        indexX=c(k,j)
        threshold_feature = (s[k] + s[k+1])/2
        ij_i = k
        ij_j = j
      }
      id2[[k1]]= rep(NA, times=(N-1))
      id2[[k1]] = idx[1:k1]
    }
    id[[j]] = id2
  }
  DD=rbind(DL,DR)
  index_left =id[[ij_j]][[ij_i]] # column, then row
  index_right = 1:N
  index_right=index_right[-index_left] 
  
  index_left=Index[sort(index_left)]
  index_right = Index[sort(index_right)]
  #IN=rbind(index_left, index_right)
  
  result=list(index_left, index_right, which_feature, threshold_feature)
  return(result)
}