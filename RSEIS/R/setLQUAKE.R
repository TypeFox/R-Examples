`setLQUAKE` <-
function(velfile="", stafile="", delfile="", fout="")
  {
    if(missing(delfile)) { delfile=NULL}
    if(missing(fout)) { fout="setup.lquake" }
    
    
    ACON = file(description =fout, open = "w")
    cat("# default setup file for 'lquake'; ps2 model and delays", file = ACON, sep = "\n")
    cat("opt = 'cps'", file = ACON, sep = "\n")
    cat("# staelv = 0", file = ACON, sep = "\n")
    cat(paste(sep=" ", "spvel = ", velfile), file = ACON, sep = "\n")

    if(is.null(delfile))
      {
        
      }
    else
      {
        cat(paste(sep=" ", "spdel = ", delfile), file = ACON, sep = "\n")
      }
    cat(paste(sep=" ", "spsta = ", stafile), file = ACON, sep = "\n")
    cat("distwt = 200.0", file = ACON, sep = "\n")
    cat("dcut = 500.0", file = ACON, sep = "\n")
    cat("amarq = 0.01", file = ACON, sep = "\n")
    cat("tol = 0.01 0.01 0.01", file = ACON, sep = "\n")
    cat("modflg = 'Y2'", file = ACON, sep = "\n")
    
    close(ACON)
    return(fout)
  }

