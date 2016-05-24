postClassProbs <-
  function(object,class=0) {
    if (!inherits(object, "randomLCA"))
      stop("Use only with 'randomLCA' objects.\n")
    if (class==0) {
      res <- as.data.frame(object$classprob)
      names(res) <- paste("Class ",1:object$nclass,sep="")
    }
    else {
      res <- data.frame(object$classprob[,class])
      names(res) <- paste("Class ",class,sep="")
    }
    res <- cbind(data.frame(object$patterns,Freq=object$freq),res)
    res
  }   