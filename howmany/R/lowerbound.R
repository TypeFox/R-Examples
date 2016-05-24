"lowerbound" <-
function(object)
  {
    if(class(object)!="howmany") stop("lowerbound requires object of class 'howmany'")
    return(object$lowerbound)
  }

