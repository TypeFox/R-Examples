"rateratio" <-
  function(x, y = NULL,
           method = c("midp", "wald"),
           conf.level = 0.95,
           rev = c("neither", "rows", "columns", "both"),
           verbose = FALSE){
    if(is.matrix(x) && !is.null(y)){stop("y argument should be NULL")}
    if(is.null(y)){
      x <- ratetable(x, rev = rev)
    } else {
      xn <- substitute(x)
      yn <- substitute(y)
      x <- ratetable(x, y, rev = rev)
      colnames(x) <- c(xn, yn)
    }
    method <- match.arg(method)
    if(method=="midp"){
      rr <- rateratio.midp(x, conf.level = conf.level,
                           verbose = verbose)
    }    
    if(method=="wald"){
      rr <- rateratio.wald(x, conf.level = conf.level,
                           verbose = verbose)
    }
    return(rr)
  }
