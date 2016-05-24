doQmapDIST <- function(x,fobj,...){
  if(!any(class(fobj)=="fitQmapDIST"))
    stop("class(fobj) should be fitQmapDIST")
  UseMethod("doQmapDIST")
}

doQmapDIST.default <- function(x,fobj,...){
### transform 'x' using quantile mapping with
### transfere function (tfun) and and parameters (par)
### specifyed in 'fobj'
### based on the parameters specifyed in 'fobj'.
###
### x : a numeric vector
### fobj : output from fitQmap (class: "fitQmap")
### ... : not used but possibly needed in further development.
  ffun <- fobj$tfun
  funpar <- fobj$par
  x <- do.call(ffun,c(list(x),as.list(funpar)))
  return(x)
}

doQmapDIST.matrix <- function(x,fobj,...){
### doQmap for matrix of time series to be
### transformed
###
### x: matrix with N columns
  if(ncol(x)!=nrow(fobj$par))
    stop("'ncol(x)' and 'nrow(fobj$par)' should be eaqual\n")  
  NN <- ncol(x)
  hind <- 1:NN
  names(hind) <- colnames(x)
  tfun <- fobj$tfun
  funpar <- fobj$par
  xx <- sapply(hind,function(i){
    tr <- try(do.call(tfun,c(list(x[,i]),as.list(funpar[i,]))),
              silent=TRUE)
    if(any(class(tr)=="try-error")){
      warning("Quantile mapping for ",names(hind)[i],
              " failed NA's produced.")
      tr <- rep(NA,nrow(x))
    }
    return(tr)
  })
  rownames(xx) <- rownames(x)
  return(xx)
}

doQmapDIST.data.frame <-  function(x,fobj,...){
### doQmap for data.frame of time series to be
### transformed
###
### x: data.frame with N columns
  x <- as.matrix(x)
  x <- doQmapDIST.matrix(x,fobj,...)
  x <- as.data.frame(x)
  return(x)
}
