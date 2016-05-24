#' @title Total, within- and between-group variances
#' 
#' @description 
#' Calculation of total, within- and between-group variance-covariance matrices
#' 
#' @param Data a numeric matrix or data frame
#' @param Group a vector of factors associated with group structure
#' @return list with the following results:
#' @return \item{Within.Var}{within-group variance-covariance matrix}
#' @return \item{Between.Var}{between-group variance-covariance matrix}
#' @return \item{Total.Var}{total variance-covariance matrix}
#' @return \item{Btween.per}{Within-group variance percentage}
#' @return \item{Btween.per}{Between-group variance percentage}
#' @seealso \code{\link{mgPCA}}, \code{\link{DGPA}}, \code{\link{DCCSWA}}, \code{\link{DSTATIS}}, \code{\link{BGC}}, \code{\link{summarize}}, \code{\link{iris}}  
#' @export
#' @references A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard (2013). General overview
#'  of methods of analysis of multi-group datasets,
#'  \emph{Revue des Nouvelles Technologies de l'Information}, 25, 108-123.
#'  
#' @examples
#' Data = iris[,-5]
#' Group = iris[,5]
#' res = TBWvariance(Data, Group)
TBWvariance <- function(Data, Group){
  
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
  if(is.null(colnames(Data))) {
    colnames(Data) = paste('V', 1:ncol(Data), sep='')
  }
  Group = as.factor(Group)
  
  
  
  rownames(Data) = Group                 #----rownames of data=groups
  M = length(levels(Group))              #----number of groups: M
  P = dim(Data)[2]                       #----number of variables: P
  n = as.vector(table(Group))            #----number of individuals in each group
  N = sum(n)                             #----number of individuals
  
  #==========================================================================
  #        	                     Outputs
  #==========================================================================
  res <- list()
  
  #=========================================================================
  #                  Within-group and between-group data
  #=========================================================================
 
  Data = scale(Data, center=TRUE, scale= FALSE) 
  dummay.matrix = model.matrix(~-1 + Group)   
  proj = dummay.matrix %*% solve(t(dummay.matrix) %*% dummay.matrix) %*% t(dummay.matrix) # proj matrix to calculate mean in each level
  res$between.data = proj %*% Data
  res$within.data  = Data - res$between.data
    
  
  res$Between.Var = var(res$between.data)
  res$Within.Var  = var(res$within.data)
  res$Total.Var   = var(Data)
  res$Within.per  = 100 * sum(diag((res$Within.Var)))  /sum(diag((res$Total.Var)))
  res$Btween.per  = 100 * sum(diag(t(res$Between.Var))) /sum(diag(t(res$Total.Var)))
  
  
  # add class
  class(res) = "TBWvariance"
  return(res)
}


#' @S3method print mgpca
print.TBWvariance <- function(x, ...)
{
  cat("\nTotal, within and between-group variances\n")
  cat(rep("-",43), sep="")
  cat("\n$Toatal.Var  ",          "Total variance")
  cat("\n$Within.Var  ",          "Within variance")
  cat("\n$Between.Var ",         "Between variance")
  cat("\n$Within.per  ",   "Within-group variance percentage")
  cat("\n$Btween.per  ",   "Between-group variance percentage")
  cat("\n")
  invisible(x)
}
