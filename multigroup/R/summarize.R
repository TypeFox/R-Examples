#' @title Summary
#' 
#' @description 
#' Summary of multigroup data in global and group parts
#' 
#' @param Data a numeric matrix or data frame
#' @param Group a vector of factors associated with group structure
#' @return list with the following results:
#' @return \item{Global.summary}{     summary of globala data}
#' @return \item{Group.summary}{      summary of group datasets}
#' @return \item{mean.between.data}{      matrix of Group mean}
#' @return \item{mean.within.data}{     matrix of group centered data}
#' @seealso  \code{\link{mgPCA}}, \code{\link{DGPA}}, \code{\link{DCCSWA}}, \code{\link{DSTATIS}}, \code{\link{BGC}}, \code{\link{TBWvariance}}, \code{\link{iris}}   
#' @export
#' @examples
#' Data = iris[,-5]
#' Group = iris[,5]
#' res = summarize(Data, Group)
summarize <- function(Data, Group){

  #=========================================================================
  #                              Checking the inputs
  #=========================================================================
  check(Data, Group)
  
  
  #=========================================================================
  #                              preparing Data
  #=========================================================================
  if (class(Data) == 'data.frame') {
    Data=as.matrix(Data)
  }
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
  
  #==========================================================================
  #      		                     Outputs
  #==========================================================================
  res <- list()
  res$Global= matrix(0, nrow=2, ncol=P)
  res$Group = list()
  #=========================================================================
  #                                 Global summary
  #=========================================================================
  
  res$Global[1,] = apply(Data, 2, mean)
  res$Global[2,] = apply(Data, 2, var)
  rownames(res$Global) = c("mean", "var") 
  colnames(res$Global) = colnames(Data) 
  
  #=========================================================================
  #                                 Group summary
  #=========================================================================
  res$Group$mean = as.matrix(aggregate(Data, by=list(Group), mean)[,-1])
  rownames(res$Group$mean) = levels(Group)
  
  res$Group$var = as.matrix(aggregate(Data, by=list(Group), var)[,-1])
  rownames(res$Group$var) = levels(Group)
  
  #=========================================================================
  #                  Mean of Within-group and between-group data
  #=========================================================================
  dummay.matrix = model.matrix(~-1 + Group)   
  proj = dummay.matrix %*% solve(t(dummay.matrix) %*% dummay.matrix) %*% t(dummay.matrix) # proj matrix to calculate mean in each level
  res$mean.between.data = proj %*% Data
  res$mean.within.data  = Data - res$mean.between.data
  
  # add class
  class(res) = "summarize"
  return(res)
}


#' @S3method print summarize
print.summarize <- function(x, ...)
{
  cat("\nsummary of multi-group data\n")
  cat(rep("-",43), sep="")
  cat("\n$Global        ",        "Global summary")
  cat("\n$Group         ",         "Group summary")
  cat("\n$between.data  ",  "Matrix of Group mean")
  cat("\n$within.data   ",   "Matrix of group centered data")
  cat("\n")
  invisible(x)
}

