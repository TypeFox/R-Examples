"EMPIRgrid" <-
function(para=NULL, deluv=0.05, verbose=FALSE, ...) {

  if(is.null(para)) {
    warning("parameters (the observed u and v) for empirical copula are NULL")
    return(NULL)
  }
  if(exists("para", para)) {
     if(length(names(para$para)) != 2) {
        warning("a data.frame having only two columns is required in the 'para$para' argument")
        return(NULL)
     }
  } else if(length(names(para)) != 2) {
     warning("a data.frame having only two columns is required in the 'para' argument")
     return(NULL)
  }

  us <- seq(0,1, by=deluv);
  vs <- seq(0,1, by=deluv);
  nu <- length(us);
  nv <- length(vs);
  cop <- matrix(nrow=nu, ncol=nv)
  if(verbose) message("Index ")
  for(i in 1:nu) {
     if(verbose) message(nu - i,"-", appendLF=FALSE)
     cop[i,] <- EMPIRcop(rep(us[i], nv), vs, para=para)
  }
  if(verbose) message("done")

  rownames(cop) <- as.character(us)
  colnames(cop) <- as.character(vs)

  zzz <- list(u=us, v=us, empcop=cop, deluv=deluv)
  return(zzz)
}


