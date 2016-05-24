#' Prediction using Random Forest or Multivariate Random Forest
#' 
#' Builds Model of Random Forest(if number of output feature is 1) or Multivariate Random Forest(if number of output
#' feature is greater than 1) using Training samples and do the prediction of testing samples using this model
#'  
#' @param trainX Input matrix of M x N, M is the number of training samples and N is the number of features 
#' @param trainY Output response of M x T, M is the number of samples and T is number of ouput Features(Response)
#' @param n_tree Number of trees in the forest
#' @param mtree Number of randomly selected features used for each split
#' @param min_leaf Minimum number of samples in the leaf node. If a node has less than equal to min_leaf samples,
#'  then there will be no splitting in that node and this node is a leaf node. 
#' @param testX Testing samples of Q x N, Q is the number of testing samples and N is the number of features(same order and
#' size used as training)
#' @return Prediction result of the Testing samples 
#' @details
#' Random Forest(RF) regression refers to ensembles of regression trees where a set of T un-pruned regression
#' trees are generated based on bootstrap sampling from the original training data. For each node, the optimal
#' feature for node splitting is selected from a random set of m features from the total N features. The selection
#' of the feature for node splitting from a random set of features decreases the correlation between different
#' trees and thus the average prediction of multiple regression trees is expected to have lower variance than
#' individual regression trees. Larger m can improve the predictive capability of individual trees but can also
#' increase the correlation between trees and void any gains from averaging multiple predictions. The bootstrap
#' resampling of the data for training each tree also increases the variation between the trees.
#' 
#' In a node with training predictor features(X) and output feature vectors(Y), node splitting is done
#' with the aim to select a feature from a random set of m features and a threshold z to partition the node 
#' into two child nodes, left node (with samples < z) and right node (with samples > z). In multivariate trees(MRF) 
#' node cost is measured as the sum of squares of the Mahalanobis 
#' distance where in univariate trees(RF) node cost is measured as the Euclidean distance.
#' 
#' After the Model of the forest is built using training predictor features(trainX) and output feature vectors(trainY),
#' the Model is used to do the prediction of features of the testing samples(testX).
#' 
#' Multivariate Random Forest(MRF) calculates prediction of all output features using one model, which is generated using
#' the training output features. While, if the output features are highly correlated then, the prediction
#' using MRF is much improved then prediction using RF. 
#' @examples
#' trainX=matrix(runif(50*100),50,100)
#' trainY=matrix(runif(50*5),50,5)
#' n_tree=5
#' mtree=10
#' min_leaf=5
#' testX=matrix(runif(10*100),10,100)
#' Command=2#2 for MRF method
#' #Prediction size is 10 x 5, where 10 is the number 
#' #of testing samples and 5 is the number of output features
#' Prediction=MultivariateRandomForest(trainX, trainY, n_tree, mtree, min_leaf, testX)
#' 
#' @references 
#' [1] Breiman, Leo. "Random forests." Machine learning 45.1 (2001): 5-32.
#' 
#' [2] Segal, Mark, and Yuanyuan Xiao. "Multivariate random forests." Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery 1.1 (2011): 80-87.
#' @export
#' 
MultivariateRandomForest <- function(trainX, trainY, n_tree, mtree, min_leaf, testX){
  theta <- function(trainX){trainX}
  results <- bootstrap::bootstrap(1:nrow(trainX),n_tree,theta) #no indics, gives number
  b=results$thetastar
  
  Variable_number=ncol(trainY)
  if (Variable_number>1){
    Command=2
  }else if(Variable_number==1){
    Command=1
  } 
  Y_HAT=matrix(  0*(1:Variable_number*nrow(testX)),  ncol=Variable_number,   nrow=nrow(testX)  )
  Y_pred=NULL
  
  for (i in 1:n_tree){
    Single_Model=NULL
    X=trainX[ b[ ,i],  ]
    Y=matrix(trainY[ b[ ,i],  ],ncol=Variable_number)
    V_inv = solve(stats::cov(Y)) # calculate the V inverse
    Single_Model=build_single_tree(X, Y, mtree, min_leaf,V_inv,Command)
    Y_pred=single_tree_prediction(Single_Model,testX,Variable_number)
    for (j in 1:Variable_number){
      Y_HAT[,j]=Y_HAT[,j]+Y_pred[,j]
    }
  }
  Y_HAT=Y_HAT/n_tree
  return(Y_HAT)
}