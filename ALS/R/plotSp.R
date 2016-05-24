plotS <- function(S, x2, out="", 
                  filename=paste("S.", out,sep=""),
                  col=vector(), cex = 1, lab="", cex.lab=1) {
  if(out=="pdf")
    pdf(filename)
  if(out=="ps")
    postscript(filename)
  par(mgp = c(2, 1, 0), mar=c(0,0,0,0), oma = c(1,0,4,0))
  if(length(lab)>0)
    par(oma=c(4,0,4,0), cex.lab=cex.lab)
  if(ncol(S) > 2)
      par(mfrow=c( ceiling( ncol(S) / 2), 2))
  if(ncol(S) == 2) # make 1 row 2 col if only plotting 2 comp.
    par(mfrow=c(2,1))
  if(ncol(S) == 1)
    par(mfrow=c(1,1))
  par(cex=cex)
  for(j in 1:ncol(S)) {
    names.side <- 0 
    if(j <= 2)
      names.side <- 3
    if( j >= ncol(S)-1 ){
      names.side <- 1 	
      par(mar=c(2,0,0,0))
    }
    if(length(col)==0)
      cl <- j
    else
      cl <- col[j]
    if(names.side != 0)	
      barplot3(S[,j], col = cl, border = cl, 
               names.arg = x2, 
               names.side = names.side, 
               names.by = 30,
               axes=FALSE)
    else 
      barplot3(S[,j], col = cl, border = cl, axes=FALSE,
               names.arg="")
  }
  if(length(lab)>0)
    mtext(lab, side = 1, outer = TRUE, line = 1, cex=cex.lab)
  if(out=="ps"||out=="pdf") dev.off()
}
