"f2fpds" <-
function(f, rate=NA) {
   if(is.na(rate)) return(f)
   if(! check.fs(f)) return(FALSE)
   if(rate <= 0) {
      warning("invalid rate ? <= 0 parameter")
      return(rep(NA,length(f)))
   }
   g <- (log(f) + rate) / rate
   g[g < 0] <- NA
   g[g > 1] <- NA
   return(g)
}

