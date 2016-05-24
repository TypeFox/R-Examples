"blomCOP" <-
function(cop=NULL, para=NULL, as.sample=FALSE, ...) {
   if(as.sample) {
      if(is.null(para)) {
         warning("Sample Blomqvist's Beta desired but para is NULL, returning NULL")
         return(NULL)
      }
      if(length(names(para)) != 2) {
        warning("para argument must be data.frame having only two columns, returning NULL")
        return(NULL)
      }
      u <- para[,1]; v <- para[,2]; n <- length(u); A <- (1+n)/2
      samBLOM <- (2/n)*(sum(as.numeric((rank(u) - A)*(rank(v) - A) >= 0))) - 1
      return(samBLOM)
   } else {
      if(is.null(cop)) {
         warning("must have copula argument specified, returning NULL")
         return(NULL)
      }
      blom <- 4*cop(0.5,0.5, para=para, ...) - 1
      return(blom)
   }
}

