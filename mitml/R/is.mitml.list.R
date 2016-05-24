is.mitml.list <- function(x){
# checks if the argument is a list of class "mitml.list"

  l <- any(class(x)=="mitml.list") & is.list(x)
  if(!l){
    return(FALSE)
  }else{
    if(any(!sapply(x,is.data.frame))) warning("Does not appear to be a list of data frames.")
    return(TRUE)
  }
}
