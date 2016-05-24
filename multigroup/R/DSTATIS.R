#' @title Dual STATIS
#' 
#' @description 
#' Dual STATIS
#' 
#' @param Data a numeric matrix or data frame
#' @param Group a vector of factors associated with group structure
#' @param ncomp number of components, if NULL number of components is equal to 2
#' @param Scale scaling variables, by defalt is False. By default data are centered within groups. 
#' @param graph should loading and component be plotted
#' @return list with the following results:
#' @return \item{Data}{original data}
#' @return \item{Con.Data}{Concatenated centered data}
#' @return \item{split.Data}{Group centered data}
#' @return \item{Group}{Group as a factor vector}
#' @return \item{RV}{The RV coefficient matrix}
#' @return \item{weights}{Vector of weights}
#' @return \item{compromise.matrix}{Compromise variance-covariance matrix}
#' @return \item{loadings.common}{Matrix of common loadings}
#' @return \item{lambda}{The specific variances of group}
#' @seealso \code{\link{mgPCA}}, \code{\link{FCPCA}}, \code{\link{DCCSWA}}, \code{\link{BGC}}, \code{\link{DGPA}}, \code{\link{summarize}}, \code{\link{TBWvariance}}, \code{\link{loadingsplot}}, \code{\link{scoreplot}}, \code{\link{iris}}  
#' @export
#' @references 
#' C. Lavit (1988). \emph{Analyse conjointe de tableaux quantitatifs}. Masson.
#'    
#' C. Lavit, Y. Escoufier, R. Sabatier and P. Traissac (1994).
#' The ACT (STATIS method). \emph{Computational Statistics & Data Analysis}, 18, 97-117. 
#' 
#' A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard (2013). General overview
#'  of methods of analysis of multi-group datasets,
#'  \emph{Revue des Nouvelles Technologies de l'Information}, 25, 108-123.
#' 
#'    
#'        
#' @examples
#' Data = iris[,-5]
#' Group = iris[,5]
#' res.DSTATIS = DSTATIS(Data, Group, graph=TRUE)
#' loadingsplot(res.DSTATIS, axes=c(1,2))
#' scoreplot(res.DSTATIS, axes=c(1,2)) 
DSTATIS <- function(Data, Group, ncomp=NULL, Scale=FALSE, graph=FALSE){
  
  
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
  Group = as.factor(Group)
  
  
  
  rownames(Data) = Group                 #---- rownames of data=groups
  M = length(levels(Group))              #----number of groups: M
  P = dim(Data)[2]                       #----number of variables: P
  n = as.vector(table(Group))            #----number of individuals in each group
  N = sum(n)                             #----number of individuals
  split.Data = split(Data,Group)         #----split Data to M parts 
  
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
  #    			                       Outputs
  #==========================================================================
  res <- list(
    Data       = Data,
    Con.Data   = Con.Data,
    split.Data = split.Data,
    Group=Group)

  #============================================================================
  #    			                      Analysis
  #============================================================================
  
  # RV coefficients
  RV = matrix(0, nrow=M, ncol=M)
  for(j in 1:M){
    for(k in 1:M){
      ab = sum( diag( cov.Group[[j]] %*% cov.Group[[k]] ) )
      aa = sum( diag( cov.Group[[j]] %*% cov.Group[[j]] ) )
      bb = sum( diag( cov.Group[[k]] %*% cov.Group[[k]] ) )
      RV[j,k] = ab/ sqrt(aa * bb)
    }
  }
  res$RV = RV

  # alpha: eigenvector of RV assocaited with the largest eigenvalue
  alpha = eigen(RV)$vectors[,1]
  alpha = abs(alpha)           
  res$weights =  alpha 

  #The compromise variance-covariance matrix Vc=sum alpha_m * V_m
  Vc = alpha[1] * cov.Group[[1]]
  for(m in 2:M){
    Vc = Vc + alpha[m] * cov.Group[[m]]
  }
  res$compromise.matrix = Vc 

 
  res$loadings.common = eigen(Vc)$vectors[, 1:ncomp]  
  rownames(res$loadings.common)=colnames(Data)
  colnames(res$loadings.common)=paste("Dim", 1:ncomp, sep="")
  
  # Lambda= A' V_m A
  res$lambda = matrix(0, nrow=M, ncol=ncomp)
  for(m in 1:M){
    res$lambda[m,] = round(diag(t(res$loadings.common) %*% cov.Group[[m]] %*% res$loadings.common),3)
  }
  rownames(res$lambda) = levels(Group)
  colnames(res$lambda) = paste("Dim", 1:ncomp, sep="")

  
  
  if(graph) {plot.mg(res)}
  
  # add class
  class(res) = c("DSTATIS", "mg")
  return(res)
}


#' @S3method print DSTATIS
print.DSTATIS <- function(x, ...)
{
  cat("\nDual STATIS\n")
  cat(rep("-",43), sep="")
  cat("\n$lambda            ", "Variance of each froup")
  cat("\n$loadings.common   ", "common loadings")
  cat("\n$Data              ", "Data set")
  cat("\n")
  invisible(x)
}

