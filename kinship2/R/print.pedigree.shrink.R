# Automatically generated from all.nw using noweb
#$Log: print.pedigree.shrink.q,v $
#Revision 1.2  2009/11/19 14:35:01  sinnwell
#add ...
#
#Revision 1.1  2009/11/17 14:39:32  sinnwell
#Initial revision
#
#Revision 1.1  2008/07/16 20:23:14  sinnwell
#Initial revision
#

print.pedigree.shrink <- function(x, ...){

  printBanner(paste("Shrink of Pedigree ", unique(x$pedObj$ped), sep=""))
 
  cat("Pedigree Size:\n")

  if(length(x$idTrimmed) > 2)
    {
      n <- c(x$pedSizeOriginal, x$pedSizeIntermed, x$pedSizeFinal)
      b <- c(x$bitSize[1], x$bitSize[2], x$bitSize[length(x$bitSize)])
      row.nms <- c("Original","Only Informative","Trimmed")
    } else {
      n <- c(x$pedSizeOriginal, x$pedSizeIntermed)
      b <- c(x$bitSize[1], x$bitSize[2])
      row.nms <- c("Original","Trimmed")
    }

  df <- data.frame(N.subj = n, Bits = b)
  rownames(df) <- row.nms
  print(df, quote=FALSE)
  
  
  if(!is.null(x$idList$unavail)) 
    cat("\n Unavailable subjects trimmed:\n", x$idList$unavail, "\n")
  
  if(!is.null(x$idList$noninform)) 
    cat("\n Non-informative subjects trimmed:\n", x$idList$noninform, "\n")
  
  if(!is.null(x$idList$affect)) 
    cat("\n Informative subjects trimmed:\n", x$idList$affect, "\n")
  
  
  ##cat("\n Pedigree after trimming:", x$bitSize, "\n")
  
  invisible()
}

