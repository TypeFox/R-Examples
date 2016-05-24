#' Split of the Parent node
#' 
#' Split of the training samples of the parent node into the child node 
#' for the feature which gives the minimum cost of splitting
#'  
#' @param X Input Training matrix of M x N, M is the number of training samples and N is the number of features
#' @param Y Output Training response of M x T, M is the number of samples and T is number of ouput Features(Response)
#' @param mtree number of randomly selected features used for each split
#' @param Index Index of training samples
#' @param V_inv Covariance matrix of Output Feature matrix
#' @param Command 1 for RF and 2 for MRF depending on the method
#' @return List with the following components:
#' \item{index_left}{Index of the samples which are in the left node after splitting}
#' \item{index_right}{Index of the samples which are in the right node after splitting}
#' \item{which_feature}{The number of the feature which gives the minimum splitting cost}
#' \item{threshold_feature}{The threshold value, which will decide a sample will go to the left node or,
#' right node. If the training sample value is less than threshold value, it will go to the left node or if greater
#' then go to the right node.}
#' @details 
#' In a node of a decision tree, there are number of samples with different features. In time of splitting, a fixed number
#' of features(mtree) has been selected randomly(that's why it is called random forest). For the splitting, node cost 
#' for all the spliting of these features are considered and whichever gives the minimum value has been selected as 
#' the splitting criteria(feature value and threshold value of the feature) of this node split. 
#' @examples
#' X=matrix(runif(20*100),20,100)
#' Y=matrix(runif(20*3),20,3)
#' mtree=5
#' Index=1:20
#' V_inv=stats::cov(Y)
#' Command=2#MRF, as number of output feature is greater than 1
#' Split_criteria=spliting(X,Y,mtree,Index,V_inv,Command) 
#' @export
spliting <- function(X,Y,mtree,Index,V_inv,Command){
  x=X[Index, ]
  if (Command==1){
    y=matrix(Y[Index, ],ncol=1)
  }else {
    y=Y[Index, ]
  }
  f = ncol(x) # number of features
  f_1 =sort(sample(f, mtree)) #randomly taken 10 features, for each splits vary
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
      D_left= Multi_D_mod(yk_left,V_inv,Command)
      D_right= Multi_D_mod(yk_right,V_inv,Command)
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