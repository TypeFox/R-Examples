"prod2COP" <-
function(u,v, cop1=NULL, para1=NULL, cop2=NULL, para2=NULL,
              para=NULL, interval=NULL, ...) {

   if(! is.null(para)) {
      cop1  <- para$cop1;    cop2  <- para$cop2
      para1 <- para$para1;   para2 <- para$para2
   }

   if(is.null(cop1)) {
        warning("must have first copula specified, returning NULL")
        return(NULL)
   }
   if(is.null(cop2)) {
        warning("must have second copula specified, returning NULL")
        return(NULL)
   }

   if(length(u) > 1 & length(v) > 1 & length(u) != length(v)) {
      warning("length u = ", length(u), " and length v = ", length(v))
      warning("longer object length is not a multiple of shorter object length, no recycling")
      return(NA)
   }
   if(length(u) == 1) {
      u <- rep(u, length(v))
   }
   else if(length(v) == 1) {
      v <- rep(v, length(u))
   }

   # d/du derCOP and d/dv derCOP2  (Nelsen, 2006, eq. 6.4.2)
   "afunc" <- function(t, U=NA, V=NA, ...) derCOP( t,V, cop=cop2, para=para2, ...) *
                                           derCOP2(U,t, cop=cop1, para=para1, ...)
   lo <- 0; hi <- 1
   if(! is.null(interval)) { lo <- interval[1]; hi <- interval[2] }
   return(sapply(1:length(u), function(i) integrate(afunc, lo,hi, U=u[i], V=v[i],
                           cop1=cop1, cop2=cop2, para1=para1, para2=para2, ...)$value))
}
