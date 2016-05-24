"giniCOP" <-
function(cop=NULL, para=NULL, by.concordance=FALSE, as.sample=FALSE, ...) {
   if(as.sample) {
      if(is.null(para)) {
         warning("Sample Gini's Gamma desired but para is NULL, returning NULL")
         return(NULL)
      }
      if(length(names(para)) != 2) {
        warning("para argument must be data.frame having only two columns, returning NULL")
        return(NULL)
      }
      R <- rank(para[,1]); S <- rank(para[,2]); n <- length(para[,1])
      samGAM <- sum(abs((n + 1 - R) - S) - abs(R - S))
      return( (2 / n^2) * samGAM)
   }

   if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }
   if(by.concordance == TRUE) {
     A <- concordCOP(cop=cop, para=para, cop2=M, ...)
     B <- concordCOP(cop=cop, para=para, cop2=W, ...)
     gini <- A + B
     return(A+B)
   }
   myC1 <- function(x, para=NULL) return(cop(x,   x, para=para, ...))
   myC2 <- function(x, para=NULL) return(cop(x, 1-x, para=para, ...))
   C1 <- integrate(myC1, 0, 1, para=para)
   C2 <- integrate(myC2, 0, 1, para=para)
   A <- C1$value; B <- C2$value
   gini <- 4*(A+B) - 2
   return(gini)
}
