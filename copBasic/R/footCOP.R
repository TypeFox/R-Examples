"footCOP" <-
function(cop=NULL, para=NULL, by.concordance=FALSE, as.sample=FALSE, ...) {
   if(as.sample) {
      if(is.null(para)) {
         warning("Sample Spearman's Footrule desired but para is NULL, returning NULL")
         return(NULL)
      }
      if(length(names(para)) != 2) {
        warning("para argument must be data.frame having only two columns, returning NULL")
        return(NULL)
      }
      samPSI <- 1 - 3*sum(abs(rank(as.numeric(para[,1])) -
                              rank(as.numeric(para[,2]))))   /   (length(para[,1])^2 - 1)
      return(samPSI)
   }

   if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }
   if(by.concordance) {
      Q <- concordCOP(cop=cop, para=para, cop2=M, ...)
      return((3/2)*Q - (1/2))
   }
   "afunc" <- function(u) { cop(u,u, para=para, ...) }
   my.int <- NULL
   try(my.int <- integrate(afunc, 0,1, ...) )
   if(is.null(my.int)) {
      warning("error encountered on numberical integration---try by.concordance instead?")
      return()
   }
   return(6*my.int$value - 2)
}
