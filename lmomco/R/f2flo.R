"f2flo" <-
function(f, pp=NA) {
   if(! check.fs(f)) return(FALSE)
   if(is.na(pp)) {
      warning("pp can not be NA")
      return(FALSE)
   } else {
      if(pp < 0 || pp > 1) {
        print("pp argument is not a valid nonexceedance probability")
        return(FALSE)
      }
   }
   zs <- (f-pp)/(1-pp)
   if(any(zs < 0) || any(zs > 1)) {
      warning("invalid nonexceedance probability after pp conditioning")
      return(FALSE)
   }
   return(zs)
}
