"uvlmoms" <- function(u,v=NULL, umv=TRUE, p=NA, type="gno", getlmoms=TRUE, ...) {
   if(is.null(v)) {
      if(length(names(u)) != 2) {
         warning("a data.frame having only two columns is required")
         return(NULL)
      } else {
         v <- u[,2]
         u <- u[,1]
      }
   }
   if(length(u) != length(v)) {
      warning("argument(s) or implied arguments u and v are unequal in length, ",
              "returning NULL")
      return(NULL)
   }

   ifelse(umv, psi <- u-v, psi <- u+v-1)
   ifelse(is.na(p), getlmoms <- TRUE, getlmoms <- FALSE)

   if(getlmoms) { # L-moments desired
     lmr <- lmomco::lmoms(psi, ...)
     lmr$ratios[2] <- NA # although numerically computable, theory
     # requires that the LCV for a distribution required to have
     # mean zero is NaN.
     return(lmr)
   }

   if(p < 0 | p >= 1/2) {
      warning("argument p is not in 0 < p < 1/2, returning NULL")
      return(NULL)
   }

   if(any(lmomco::dist.list() == type)) {
      # Fitting a univariate distribution by method of L-moments
   	  lmr  <- lmomco::lmoms(psi, ...)
   	  para <- lmomco::lmom2par(lmr, type)
         A <- lmomco::qlmomco(1-p, para)
         B <- lmomco::qlmomco(  p, para)
         C <- lmomco::qlmomco(1/2, para)
      skewness <- A - 2*C + B / (A - B)
      return(skewness)
   } else {
   	  # Using the empirical distribution function
   	  type <- as.integer(type)
      if(is.na(type) | type < 1 | type > 9) {
         warning("argument 'type' as an integer is not in [1, 9] as quantile() will ",
                 "require, returning NULL")
         return(NULL)
      }
      A <- quantile(psi, probs=1-p, names=FALSE, type=type)
      B <- quantile(psi, probs=  p, names=FALSE, type=type)
      C <- quantile(psi, probs=1/2, names=FALSE, type=type)
      skewness <- A - 2*C + B / (A - B)
      return(skewness)
   }
}


"uvskew" <- function(u,v=NULL, umv=TRUE, p=0.05, type=6, getlmoms=FALSE, ...) {
   return(uvlmoms(u,v=v, umv=umv, p=p, type=type, getlmoms=FALSE, ...))
}



