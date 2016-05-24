.convert.result.format<-function(x, resname="res"){
  if(is.list(x))
     {rdim <- length(x)
      cnames <- names(x[[1]])
      if(is.null(cnames))
          cnames <- paste(abbreviate(resname),1:cdim,sep=".")
      else
          cnames <- paste(abbreviate(resname),abbreviate(names(x[[1]])),sep=".")
      x <- data.frame(matrix(unlist(x),nrow=rdim, byrow=TRUE))
      colnames(x) <- cnames
     }
  else if(is.matrix(x))
     {x <- t(x)
      cdim <- ncol(x)
      cnames <- colnames(x)
      if(is.null(cnames))
          cnames <- paste(abbreviate(resname),1:cdim,sep=".")
      else
          cnames <- paste(abbreviate(resname),abbreviate(names(x[[1]])),sep=".")
      x <- data.frame(x)
      colnames(x) <- cnames
     }
  else if (!is.data.frame(x))
     {x <- data.frame(x)
      names(x) <- abbreviate(resname)
     }
  return(x)
}


