"kfuncCOPlmoms" <-
function(cop=NULL, para=NULL, nmom=5, begin.mom=1, ...) {
   if(begin.mom > nmom) {
       warning("begin.mom can not be greater than number of L-moments ",
               "requested by nmom, resetting to nmom")
       begin.mom <- nmom
   }

   if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }

   L <- R <- vector(mode="numeric", length=nmom)
   L[1:nmom] <- R[1:nmom] <- NA
   for(r in begin.mom:nmom) L[r] <- kfuncCOPlmom(r, cop=cop, para=para, ...)
   R[2] <- L[2]/L[1]
   begr <- ifelse(begin.mom > 3, begin.mom, 3)
   R[begr:nmom] <- sapply(begr:nmom,
              function(r) { ifelse(is.na(L[2]), return(NA), return(L[r]/L[2]))})
   z <- list(lambdas=L, ratios=R,
             source="kfuncCOPlmoms")
   return(z)
}


"kfuncCOPlmom" <-
function(r, cop=NULL, para=NULL, ...) {
   if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }
   if(r < 1) { warning("r < 1, returning NA"); return(NA) }
   if(r == 1) {
      tmpA <- NULL
      try(tmpA <- integrate(function(t) { 1-kfuncCOP(t,cop=cop,para=para,...) }, 0,1))
      if(is.null(tmpA)) {
         warning("could not integrate for the mean of the Kendall Function")
         return(NA)
      }
      return(tmpA$value)
   }
   # Technically, R appears to operate correctly to get the mean even when sapply'ing
   # the 0:(r-2) for r == 1. However, it is about 100 percent slower then just integrating
   # the above condition for r == 1.
   "sfunc" <- function(j) {
      "afunc" <- function(z, j) { Fz <- kfuncCOP(z, cop=cop, para=para, ...)
                                  return( Fz^(r-j-1) * (1-Fz)^(j+1) ) }
      tmpA <- NULL
      try(tmpA <- integrate(afunc, 0, 1, j=j))
      if(is.null(tmpA)) {
         warning("could not integrate for the mean of the Kendall Function")
         return(NA)
      }
      tmpA <- tmpA$value
      tmpB <- (-1)^j * exp(lchoose(r-2,j) + lchoose(r, j+1))
      return(tmpA*tmpB)
   }
   tmp <- sum(sapply(0:(r-2), sfunc)) / r
   return(tmp)
}

