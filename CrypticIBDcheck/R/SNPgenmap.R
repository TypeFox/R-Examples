SNPgenmap <- function(physmap, chromosomes) {
   genmap <- rep(NA, length(chromosomes))
   for (chr in 1:22) {
      physpos <- RutgersMapB36[[paste("chr", chr, sep="")]]$Build36_map_physical_position
      genmappos<-RutgersMapB36[[paste("chr", chr, sep="")]]$Sex.averaged_map_position
      chrmap <- approxfun(physpos,genmappos)
      ind <- which(chromosomes == chr)
      genmap[ind] <- chrmap(physmap[ind])
   }
   return(genmap)
}


