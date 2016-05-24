progRep <- function(info,verbInt,tt,tmax) {
   eps <- .Machine$double.eps
   st.now <- Sys.time()
   elt <- st.now - info$st.last
   units(elt) <- "secs"
   elt <- as.vector(elt)
   if(elt > verbInt-eps) {
      pct <- max(0,round(100*(tmax - tt)/tmax,1))
      pct <- sprintf("%4.1f",pct)
      cat("Approximately",pct,"percent of [0,tmax] left to cover.\n")
      telt <- st.now - info$st.first
      overmin <- telt > as.difftime(1,units="mins")
      if(overmin) {
         units(telt) <- "mins"
      } else {
         units(telt) <- "secs"
      }
      telt <- sprintf("%5.2f",round(as.vector(telt),2))
      cat("Total *elapsed* (wall) time approximately ",telt,
           if(overmin) " minutes.\n" else " seconds.\n",sep="")
      info$st.last <- st.now
   }
   invisible()
}
