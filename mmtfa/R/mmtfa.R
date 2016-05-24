mmtfa <- function(x, Gs=1:4, Qs=1:2, clas=0, init="kmeans", scale=TRUE, models = "all", 
                  dfstart=50, dfupdate="approx", gauss=FALSE, eps=0.05, known=NULL, parallel.cores=FALSE){
  modgen <- modelgen()
  modold  <- modgen$modold
  p <- ncol(x)
  n <- nrow(x)
  if(any(Qs>round(p/2))){
    Qs <- Qs[Qs<=round(p/2)]
    if(length(Qs)<1){
      stop("No value of Qs meets criteria: Number of factors cannot exceed half dimension size.")
    }
    else{
      warning(paste("Number of factors cannot exceed half dimension size: Q =", paste(Qs[Qs<=round(p/2)], collapse=", "), "run"))
    }
  }
  
  if(gauss){
    dfstart <- 200
    dfupdate <- FALSE
    models <- modgen$dfconstrained
  }
  if(!(dfupdate %in% c(FALSE, "approx", "numeric"))){
    warning("dfupdate argument not recognized, degrees of freedom not updated")
    dfupdate <- FALSE
  }
  if(scale){
    x <- scale(x, center=TRUE, scale=TRUE)
  }
  
  if(length(models)==1){
    if(models=="all"){
      models <- modgen$allmodels
    }
    else{
      if(models=="dfconstrained"){
        models <- modgen$dfconstrained 
      }
      else{
        if(models=="dfunconstrained"){
          models <- modgen$dfunconstrained 
        }
        else{
          if(!all(models %in% modgen$allmodels)){
            stop(paste("Unknown model(s) specified:", models[!(models %in% modgen$allmodels)]))
          }
        }
      }
    }
  }
  else{
    if(!all(models %in% modgen$allmodels)){
      stop(paste("Unknown model(s) specified:", models[!(models %in% modgen$allmodels)]))
    }
  }
  
  
  if(is.logical(parallel.cores)){
    if(parallel.cores){
      modrun <- mmtfa.parallel(x, Gs, Qs, clas, init, scale, models, 
                               dfstart, dfupdate, gauss, eps, known, numcores=NULL)
    }
    else{
      modrun <- mmtfaEM(x, Gs, Qs, clas, init, scale, models, dfstart, dfupdate, gauss, eps, known)
    }
  }
  else{ modrun <- mmtfa.parallel(x, Gs, Qs, clas, init, scale, models, 
                                 dfstart, dfupdate, gauss, eps, known, numcores=parallel.cores) }
  
  return(modrun)
}