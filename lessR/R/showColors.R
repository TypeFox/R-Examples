showColors <- 
function(file="colors.pdf", color=NULL) {

  if (!is.null(color)) file <- paste("colors.", color, ".pdf", sep="")
  pdf(file=file)

  old.par <- par(no.readonly = TRUE)  # save all modifiable settings
  on.exit(par(old.par))
  par(mfrow=c(5,6), mgp=c(0,1,0))
  
  if (is.null(color))
    clr <- colors()
  else
    clr <- grep(color, colors(), value = TRUE)
    
  h <- 1
  for (i in 1:length(clr))
    barplot(h, col=clr[i], main=clr[i], sub=toString(col2rgb(clr[i])), 
       cex.main=.95, axes=FALSE, border=NA)
  
  .showfile(file, "file  of R colors")
  
  dev.off()

}
