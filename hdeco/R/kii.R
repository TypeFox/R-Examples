"kii" <-
function(BE1,BE2,fnev="aentrop.txt",sisi=3) {
  BE2[BE1==0] <- 0
  KI <- sum(BE1*BE2)
  IK <- length(BE1)
  cat(c("Decomp:",signif(KI,sisi)),file=fnev,append=T)
  JEL <- "  = "
  for (ik in 1:IK) {
    cat(c(JEL,signif(BE1[ik],sisi),"*",signif(BE2[ik],sisi)),file=fnev,append=T,sep="")
    JEL <- " + "
  }
  cat(" ",fill=T,file=fnev,append=T)
}

