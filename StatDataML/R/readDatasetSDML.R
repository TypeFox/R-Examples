readDatasetSDML <- function(x)
{
  if(is.null(x)) return(NULL)
  
  ## a dataset contains either a list or an array 
  if(!is.null(x[["list"]]))
    return(readListSDML(x[["list"]]))
    
  if(!is.null(x[["array"]]))
    return(readArraySDML(x[["array"]]))
}
