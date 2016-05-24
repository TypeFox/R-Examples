model.matrix.msc <- function (object, newdata, ...) 
{
    #try to extract model data from newdata
    if(missing(newdata) || is.null(newdata)){
      x <- object$x
    }
    else{
      nr <- nrow(newdata)
      if(is.null(nr)){
        newdata <- matrix(newdata, nrow=1)
      }
      if(ncol(newdata) == ncol(object$x)){
        x <- as.matrix(newdata)
      } 
      else{
        x = matrix(nrow = nrow(newdata), ncol = ncol(object$x))
        cn <- colnames(object$x)
        cnd <- colnames(newdata)
        for(i in 1: ncol(object$x)){
          tmp <- newdata[cnd == cn[[i]] ]
          if(is.null(tmp)){
           stop("newdata does not match model")
          }
          x[,i] = tmp[[1]]
        }
      }
    }
    colnames(x) <- colnames(object$x)
    x
}
