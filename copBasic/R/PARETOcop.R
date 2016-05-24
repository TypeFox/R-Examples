"PARETOcop" <- function(u,v, para=NULL, ...) {
   if(para < 0) {
      warning("parameter must be 0 <= Theta < Infinity [well numerically 100]")
      return(NA)
   }
   if(length(u) > 1 & length(v) > 1 & length(u) != length(v)) {
      warning("length u = ",length(u), " and length v = ",length(v))
      warning("longer object length is not a multiple of shorter object length, ",
              "no recycling")
      return(NA)
   }
   # The extra hassle of vectorization made here is to handle situations
   # in which nested integrals are used where uneven vectors can be passed
   if(length(u) == 1) {
      u <- rep(u, length(v))
   } else if(length(v) == 1) {
      v <- rep(v, length(u))
   }

   if(para == 0)   return(P(u,v))
   if(para < 1E-8) return(P(u,v)) # from rhoCOP 
   if(para > 100)  return(M(u,v)) # from rhoCOP 

   return( ((1-u)^(-para) + (1-v)^(-para) - 1)^(-1/para) + u + v - 1 )
}
