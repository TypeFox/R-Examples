predict.PCAmix <- function(object,X.quanti=NULL,X.quali=NULL,...)
{
  pca<-object
  if (!inherits(pca, "PCAmix")) 
    stop("use only with \"PCAmix\" objects")
  
  rec <- recod(X.quanti,X.quali)
  Y <- rec$Y
  n <- rec$n
  beta <- pca$coef
  
  ncol_tot<-length(beta[[1]])-1
  #test if the subsample of individual to predict has all the levels of the categ var
  if (ncol_tot!=ncol(Y)){
    Ymodif<-matrix(0,nrow=n,ncol=ncol_tot)
    colname.part<-colnames(Y)
    colname.tot<-names(beta[[1]][-1,])
    colnames(Ymodif)<-colname.tot
    Ymodif[,colname.part]<-Y[,colname.part]
    Y<-Ymodif  
  }
    
  if (!is.null(X.quanti)) 
    {
    label <- rownames(X.quanti)
    n1 <- nrow(X.quanti)
    p1 <- ncol(X.quanti)
    if (p1 != pca$rec$p1) stop("The number of variables in X.quanti must be the same than in the learning set")
  }
  if (!is.null(X.quali))
    {
    label <- rownames(X.quali)
    n2 <- nrow(X.quali)
    p2 <- ncol(X.quali)
    if (p2 != pca$rec$p2) stop("The number of variables in X.quali must be the same than in the learning set")
  }
  if (!is.null(X.quanti)&& !is.null(X.quali))
    {
    if (n1 != n2) stop("The number of objects in X.quanti and X.quali must be the same")
    if (sum(rownames(X.quali)!=rownames(X.quanti))!=0) stop("The names of the objects in X.quanti and X.quali must be the same")
  }
  
  coord <- matrix(,n,length(beta))
  for (g in 1: length(beta)) coord[,g] <-Y %*% beta[[g]][-1] +  beta[[g]][1]
  
  if (colnames(pca$sqload)[1]=="dim1.rot")  
    colnames(coord) <- paste("dim", 1:length(beta), sep = "",".rot")
  else
    colnames(coord) <- paste("dim", 1:length(beta), sep = "")
  rownames(coord) <- label
  return(coord)			
}
