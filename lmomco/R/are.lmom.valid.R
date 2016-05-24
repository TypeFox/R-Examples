"are.lmom.valid" <-
function(lmom) {

   if(length(lmom$L1) == 0) { # convert to named L-moments
     lmom <- lmorph(lmom)     # nondestructive conversion!
   }

   # The early return trues are for situations in which the higher moments
   # are simply not available--say from computing the l-moments of a distribution
   if(is.null(lmom$L2) || is.na(lmom$L2))    return(TRUE)
   if(lmom$L2 <= 0)        return(FALSE)
   if(is.null(lmom$TAU3) || is.na(lmom$TAU3))  return(TRUE)
   if(abs(lmom$TAU3) >= 1) return(FALSE)
   if(is.null(lmom$TAU4) || is.na(lmom$TAU4))  return(TRUE)
   if(lmom$TAU4 < (0.25*(5*lmom$TAU3^2 - 1)) || lmom$TAU4 >= 1) return(FALSE)
   if(is.null(lmom$TAU5) || is.na(lmom$TAU5))  return(TRUE)
   if(abs(lmom$TAU5) >= 1) return(FALSE)
   return(TRUE)
} 

