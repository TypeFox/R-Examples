#' @title  Transter function
#' 
#' @description 
#'  Transfer data to M groups of individuals
#' 
#' @param Data a numeric (quantitative) matrix or data frame
#' @param Group a vector of factors associated with group structure
#' @return list with the following results:
#' @return \item{Data.group}{original data}
#' @seealso \code{\link{split}}, \code{\link{iris}} 
#' @keywords internal
transfer_Group = function(Data, Group){
   
  if (class(Data) == 'data.frame') {
    Data = as.matrix(Data)
  }
  rownames(Data) = Group                 #---- rownames of data=groups
  M = length(levels(Group))
  P = dim(Data)[2]
  n = as.vector(table(Group))
  N = sum(n)
  
  #============================================================================
  #  	       	                          output
  #============================================================================
  Data.group=list()
  
  #============================================================================
  #					                          split data
  #============================================================================
  Xm = split(Data, Group)  #Split Data
  grouprname = split(rownames(Data), Group)
  
  for(m in 1:M){  # dataset in form of matrix
    
    if(class(Xm[[m]])!="data.frame"){
      Xm[[m]] = matrix(Xm[[m]], ncol=P)}
    colnames(Xm[[m]]) = colnames(Data)
    rownames(Xm[[m]]) = grouprname[[m]]
  }
  
  
  Data.group=Xm
  
  
  return(Data.group)
}
