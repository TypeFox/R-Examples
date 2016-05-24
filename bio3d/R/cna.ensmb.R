# ensemble CNA calculation optimized with multicore
cna.ensmb <- function(cij, ..., ncore = NULL) {
   ncore <- setup.ncore(ncore)

   cijs <- cij 
   if("all.dccm" %in% names(cijs)) cijs <- cijs$all.dccm 
   if(is.array(cijs) && length(dim(cijs))==3)
      cijs <- do.call("c", apply(cijs, 3, list))
   if(is.list(cijs)) {
      net <- mclapply(cijs, cna.dccm, ...)
   } else {
      warning("cijs should be matrix, array(dim=3), or list")
      net <- NULL 
   }
   return(net)
}
