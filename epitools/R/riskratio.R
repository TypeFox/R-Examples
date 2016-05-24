"riskratio" <-
  function(x, y = NULL,
           method = c("wald", "small", "boot"),
           conf.level = 0.95,
           rev = c("neither", "rows", "columns", "both"),
           correction = FALSE,
           verbose = FALSE,
           replicates = 5000){
    if(is.matrix(x) && !is.null(y)){stop("y argument should be NULL")}
    if(is.null(y)){
      x <- epitable(x, rev = rev)
    } else {
      x <- epitable(x, y, rev = rev)
    }
    method <- match.arg(method)
    if(method=="wald"){
      rr <- riskratio.wald(x, conf.level = conf.level,
                           verbose = verbose, correction = correction)
    }
    if(method=="small"){
      rr <- riskratio.small(x, conf.level = conf.level,
                            verbose = verbose, correction = correction)
    }
    if(method=="boot"){
      rr <- riskratio.boot(x, conf.level = conf.level,
                           replicates = replicates,
                           verbose = verbose, correction = correction)
    }
    rr
  }
