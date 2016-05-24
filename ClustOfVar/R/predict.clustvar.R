#' @S3method predict clustvar
#' @export
#' @name predict.clustvar
#' @title Scores of new objects on the synthetic variables of a given partition
#' @description A partition of variables obtained with kmeansvar or with cutreevar is given in input.
#' Each cluster of this partition is associated with a synthetic variable which is a linear combination of the variables of the cluster.
#' The coefficients of these k linear combinations (one for each cluster) are used here to calculate new scores of a objects described in a new dataset (with the same variables).
#' The output is the matrix of the scores of these new objects on the k synthetic variables.
#' @param part  an object of class clustvar
#' @param X.quanti  numeric matrix of data for the new objects
#' @param X.quali  a categorical matrix of data for the new objects
#' @return Returns the matrix of the scores of the new objects on the k syntetic variables of the k-clusters partition given in input.
#' @author Marie Chavent \email{marie.chavent@@math.u-bordeaux1.fr}, Vanessa Kuentz, Benoit Liquet, Jerome Saracco
#' @examples  
#' data(wine)
#' n <- nrow(wine)
#' sub <- sample(1:n,10)
#' X.quanti <- wine[sub,c(3:29)] #learning sample
#' X.quali <- wine[sub,c(1,2)] 
#' part <-kmeansvar(X.quanti,X.quali,init=5)
#' X.quanti.t <- wine[-sub,c(3:29)] 
#' X.quali.t <- wine[-sub,c(1,2)] 
#' new <- predict(part,X.quanti.t,X.quali.t)


predict.clustvar <- function(object,X.quanti=NULL,X.quali=NULL,...)
{
  part <- object
  if (!inherits(part, "clustvar")) 
    stop("use only with \"clustvar\" objects")
  
  pfin <- part$cluster
  indexj <- part$rec$indexj
  indexg<-NULL
  for (i in 1:length(indexj))
    indexg[i]<-pfin[indexj[i]] #ds le cas quanti, indexg=pfin
  if (!is.null(X.quali)) 
    G <- recodqual(X.quali) else G <- NULL
  if (!is.null(X.quanti)) 
    Y1 <- as.matrix(X.quanti) else Y1 <- NULL
  Y <- cbind(Y1,G)
  n <- nrow(Y)
  beta <- part$coef
  
  if (ncol(Y)!=ncol(part$rec$Y))
    stop("The number of categories in the learning set is different than in X.quali")
  
  if (!is.null(X.quanti)) 
  {
    label <- rownames(X.quanti)
    n1 <- nrow(X.quanti)
    p1 <- ncol(X.quanti)
    if (p1 != part$rec$p1) stop("The number of variables in X.quanti must be the same than in the learning set")
  }
  if (!is.null(X.quali))
  {
    label <- rownames(X.quali)
    n2 <- nrow(as.matrix(X.quali))
    p2 <- ncol(as.matrix(X.quali))
    if (p2 != part$rec$p2) stop("The number of variables in X.quali must be the same than in the learning set")
  }
  if (!is.null(X.quanti)&& !is.null(X.quali))
  {
    if (n1 != n2) stop("The number of objects in X.quanti and X.quali must be the same")
    if (sum(rownames(X.quali)!=rownames(X.quanti))!=0) stop("The names of the objects in X.quanti and X.quali must be the same")
  }
  
 
  scores <- matrix(,nrow(Y),length(beta))
  for (g in 1: length(beta))
  {
    Yg <- as.matrix(Y[,which(indexg==g)])
    scores[,g] <-Yg %*% beta[[g]][-1] +  beta[[g]][1]
  }
  colnames(scores) <- paste("cluster", 1:length(beta), sep = "")
  rownames(scores) <- label
  return(scores)
}