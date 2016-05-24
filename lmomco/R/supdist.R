"supdist" <- function(para, trapNaN=FALSE, delexp=0.5, paracheck=TRUE, ...) {
   if(paracheck) {
         if(! are.par.valid(para, ...)) {
         warning("parameter object seems to to be invalid, returning NULL")
         return()
      }
   }
   flo <- flo.org <- 0
   lo <- par2qua(flo.org, para, paracheck=FALSE, ...)
   lo.e <- .Machine$sizeof.longdouble
   if(is.na(lo)) {
      trapNaN <- TRUE
      lo <- NaN
   }
   if(trapNaN) {
      if(is.nan(lo)) {
         while(is.nan(lo)) {
            if(lo.e == 4) break
            flo <- flo.org + 10^(-lo.e)
            lo <- par2qua(flo, para, paracheck=FALSE, ...)
            if(is.na(lo)) lo <- NaN
            if(! is.nan(lo)) break
            lo.e <- lo.e - delexp
         }
      } else {
         lo.e <- NA
      }
   } else {
      lo.e <- NA
   }

   fhi <- fhi.org <- 1
   hi <- par2qua(fhi.org, para, paracheck=FALSE, ...)
   hi.e <- .Machine$sizeof.longdouble
   if(is.na(hi)) {
      trapNaN <- TRUE
      hi <- NaN
   }
   if(trapNaN) {
      if(is.nan(hi)) {
         while(is.nan(hi)) {
            if(hi.e == 3) break
            fhi <- fhi.org - 10^(-hi.e)
            hi <- par2qua(fhi, para, paracheck=FALSE, ...)
            if(is.na(hi)) hi <- NaN
            if(! is.nan(hi)) break
            hi.e <- hi.e - delexp
         }
      } else {
         hi.e <- NA
      }
   } else {
      hi.e <- NA
   }
   if(lo > hi) {
      warning("SERIOUS FAILURE IN supdist(), results are unreliable")
   }
   zz <- list(type       = para$type,
              support    = c(          lo,            hi),
              nonexceeds = c(         flo,           fhi),
              fexpons    = c(         -lo.e,         -hi.e),
              finite     = c(is.finite(lo), is.finite(hi)),
              source     = "supdist")
   return(zz)
}
