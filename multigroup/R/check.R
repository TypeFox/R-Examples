check <- function(Data, Group, ncomp, Scale) {
  
 
  Group = as.factor(Group)

  if(!is.data.frame(Data) & !is.matrix(Data))
    stop("\nOops the class of Data must be data.frame or matrix")
  
  # if(has_missing(Data))
  if(sum(is.na(Data)) != 0)
    stop("Oops there are missing values")
  
  
  if(nrow(Data) != length(Group))
    stop("\nOops the length of group is different of number of rows in data")
  
  
  if(length(levels(Group)) == 1)
    stop("\nGroups must have more than one level") 
  
  TRUE
}
