"fpds2f" <-
function(fpds, rate=NA) {
   if(is.na(rate)) return(fpds)
   fpds.check <- fpds[! is.na(fpds)]
   if(! check.fs(fpds.check)) return(FALSE)
   if(rate <= 0) {
       warning("invalid rate ? <= 0 parameter")
       return(rep(NA,length(fpds)))
   }
   g <- exp(-1*rate*(1 - fpds))
   g[g < 0] <- NA
   g[g > 1] <- NA
   return(g)
}
