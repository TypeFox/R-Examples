paths <-
function(object, ...)
{
  if(class(object)=="gets" || class(object)=="isat"){
    return(object$paths)
  }else{
    cat("The object does not belong to the 'gets' nor 'isat' class\n")
  }
}
