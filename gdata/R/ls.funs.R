ls.funs <- function (...)
  {
    mycall <- match.call()
    mycall[[1]] <- as.name("ls")
    nameList <- eval.parent(mycall)
    if(length(nameList)>0)
      {
        funcFlags <- sapply( nameList, function(x) is.function(get(x)) )
        return(nameList[funcFlags])
      }
    else
      return( list() )
  }
                        
