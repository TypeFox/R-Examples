# Automatically generated from all.nw using noweb
#$Log: plot.pedigree.shrink.q,v $
#Revision 1.4  2010/09/03 21:12:16  sinnwell
#use shrunk "avail" vector for the colored labels
#
#Revision 1.3  2009/11/19 14:57:18  sinnwell
#*** empty log message ***
#
#Revision 1.2  2009/11/17 23:09:51  sinnwell
#updated for ped object
#
#Revision 1.1  2008/07/16 20:23:38  sinnwell
#Initial revision
#

plot.pedigree.shrink <- function(x, bigped=FALSE, title="", 
                                 xlegend="topright", ...){

  ##  Plot pedigrees, coloring subjects according
  ##   to availability, shaded by affected status used in shrink

  if (bigped == FALSE) {
    tmp <- plot(x$pedObj, col = x$avail + 1,keep.par=T)
  }
  else {
    tmp <- plot.pedigree(x$pedObj, align = FALSE, packed = FALSE, 
                col = x$avail + 1, cex = 0.5, symbolsize = 0.5,keep.par=T)
  }
  
  legend(x = xlegend, legend = c("DNA Available", "UnAvailable"), 
         pch = c(1, 1), col = c(2, 1), bty = "n", cex=.5)
  title(paste(title, "\nbits = ", x$bitSize[length(x$bitSize)]),cex.main=.9)
  invisible(tmp)
} 

