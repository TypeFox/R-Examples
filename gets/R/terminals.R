terminals <-
function(object, ...)
{
  if(class(object)=="gets" || class(object)=="isat"){
    return(object$terminals)
  }else{
    cat("The object does not belong to the 'gets' nor 'isat' class\n")
  }
}
