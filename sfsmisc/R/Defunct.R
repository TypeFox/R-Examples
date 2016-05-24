### Functions moved from  ./Deprecated.R
###                       ~~~~~~~~~~~~~~~
###--- remove things from here to ../Old_Defunct/ex-Deprecated.R
###      ====                  == ==============================

pl.ds <- function(...) {
    stop("pl.ds() has been renamed to  plotDS() and is defunct now.\n",
          "Please change your code to use the new name")
  plotDS(...)
}

p.pllines <- function(x,y,group,lty=c(1,3,2,4),...)
{
  ## Purpose:   lines according to group
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 21 Jun 93, 15:45

  stop("p.pllines() is defunct: in R, use",
          "plot(x,y, lty=group, type='l', ...)")

  plot(x,y,type="n",...)
  ngr <- max(group)
  for (gg in 1:ngr) {
    ii <- group==gg & !is.na(x) & !is.na(y)
    if(sum(ii)) lines(x[ii],y[ii],lty=lty[1+(gg-1)%%length(lty)])
  }
}
