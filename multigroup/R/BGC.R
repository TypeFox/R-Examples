#' @title Between Group Comparison
#' 
#' @description 
#' Between Group Comparison (BGC)
#' 
#' @param Data a numeric matrix or data frame
#' @param Group a vector of factors associated with group structure
#' @param numc number of components assocaited with PCA on each group
#' @param ncomp number of components, if NULL number of components is equal to 2
#' @param Scale scaling variables, by defalt is FALSE. By default data are centered within groups
#' @param graph should loading and component be plotted
#' @return list with the following results:
#' @return \item{Data       }{Original data}
#' @return \item{Con.Data       }{Concatenated centered data}
#' @return \item{split.Data       }{Group centered data}
#' @return \item{Group       }{Group as a factor vector}
#' @return \item{loadings.common       }{Matrix of common loadings}
#' @return \item{lambda       }{The specific variances of groups}
#' @return \item{exp.var       }{Percentages of total variance recovered associated with each dimension }
#' @seealso \code{\link{mgPCA}}, \code{\link{FCPCA}}, \code{\link{DCCSWA}}, \code{\link{DSTATIS}}, \code{\link{DGPA}}, \code{\link{summarize}}, \code{\link{TBWvariance}}, \code{\link{loadingsplot}}, \code{\link{scoreplot}}, \code{\link{iris}}  
#' @export
#' 
#' 
#' @references W. J. Krzanowski (1979). Between-groups comparison of principal components,
#'  \emph{Journal of the American Statistical Association}, 74, 703-707.
#'  
#' @references A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard (2013). General overview
#'  of methods of analysis of multi-group datasets,
#'  \emph{Revue des Nouvelles Technologies de l'Information}, 25, 108-123.
#' 
#'  
#' @references A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard (2013). Analyses
#' factorielles de donnees structurees en groupes d'individus,
#' \emph{Journal de la Societe Francaise de Statistique}, 154(3), 44-57.
#' 
#'    
#'     
#'      
#'        
#' @examples
#' Data = iris[,-5]
#' Group = iris[,5]
#' res.BGC = BGC(Data, Group, graph=TRUE)
#' loadingsplot(res.BGC, axes=c(1,2))
#' scoreplot(res.BGC, axes=c(1,2)) 
BGC <- function(Data, Group, numc=NULL, ncomp=NULL, Scale=FALSE, graph=FALSE){
  
  #=========================================================================
  #                             1. Checking the inputs
  #=========================================================================
  check(Data, Group)
  
  
  #=========================================================================
  #                              2. preparing Data
  #=========================================================================
  if (class(Data) == 'data.frame') {
    Data=as.matrix(Data)
  }
  if(is.null(ncomp)) {ncomp=2}  
  if(is.null(colnames(Data))) {
    colnames(Data) = paste('V', 1:ncol(Data), sep='')
  }
  
  if(is.null(ncomp)) {ncomp=2}  
  Group = as.factor(Group)
  
  
  
  rownames(Data) = Group                 #---- rownames of data=groups
  M = length(levels(Group))              #----number of groups: M
  P = dim(Data)[2]                       #----number of variables: P
  n = as.vector(table(Group))            #----number of individuals in each group
  N = sum(n)                             #----number of individuals
  split.Data = split(Data,Group)         #----split Data to M parts 
  if(is.null(numc)) {numc=min(M-1, P)}   
  
  
  # centering and scaling if TRUE
  for(m in 1:M){  
    split.Data[[m]] = matrix(split.Data[[m]], nrow=n[m])
    split.Data[[m]] = scale(split.Data[[m]], center=TRUE, scale=Scale)
  }

  # concatinated dataset by row as groups
  Con.Data = split.Data[[1]]  
  for(m in 2:M) {
    Con.Data = rbind(Con.Data, split.Data[[m]])
  }
  rownames(Con.Data) = Group
  colnames(Con.Data) = colnames(Data)
  

  # Variance-covariance matrix for each group
  cov.Group = vector("list", M)
  for(m in 1:M){    
    cov.Group[[m]] = t(split.Data[[m]]) %*% split.Data[[m]] / n[m]
  }
  
  #==========================================================================
  #    			                     Outputs
  #==========================================================================
  res <- list(
    Data       = Data,
    Con.Data   = Con.Data,
    split.Data = split.Data,
    Group=Group)
  
  #==========================================================================
  #      		                       Method
  #==========================================================================
  #----------------- Singular Value Decomposition of a Matrix: X = U D L', 
  # selected number of components rk
  L=vector("list",M)
  for(m in 1:M){
    SVD = svd(split.Data[[m]])
    L[[m]] = SVD$v[,1:numc]
  }
  
  #----------------- H matrix:  H=sum (L'*L)
  H = matrix(0,P,P)
  for(m in 1:M){
    H = H + L[[m]] %*% t(L[[m]])
  }
  #------------- common loadings
  W = eigen(H)$vectors[,1:ncomp]  
  res$loadings.common=  W   
  rownames(res$loadings.common) = colnames(Data)
  colnames(res$loadings.common) = paste("Dim", 1:ncomp, sep="")
  
  
  #---------------- variance of each loading: lambda = t(common loading)*(t(Xm)* Xm) * common loading
  # variance of each loading # lambda = t(common loading)*(t(Xm)* Xm) * common loading
  lambda = matrix(0, nrow=M, ncol=ncomp)
  for(m in 1:M){
    lambda[m,] = round(diag(t(W) %*% cov.Group[[m]] %*% W),3)
  }
  res$lambda = lambda
  rownames(res$lambda) = levels(Group)
  colnames(res$lambda) = paste("Dim", 1:ncomp, sep="")
  
  
  #
  exp.var = matrix(0,M,ncomp)
  for(m in 1:M){
    exp.var[m,] = 100 * lambda[m,]/ sum(diag(cov.Group[[m]]))
  }
  res$exp.var = exp.var
  rownames(res$exp.var) = levels(Group)
  colnames(res$exp.var) = paste("Dim", 1:ncomp, sep="")
  

  #============================================================================
  if(graph) {plot.mg(res)}
  
  # add class
  class(res) = c("BGC", "mg")
  return(res)
}


#' @S3method print BGC
print.BGC <- function(x, ...)
{
  cat("\nBetween Group Comparison\n")
  cat(rep("-",43), sep="")
  cat("\n$loadings.common   ", "common loadings")
  cat("\n$Data              ", "Data set")
  cat("\n")
  invisible(x)
 }

  