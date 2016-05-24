
##==============================================================================
# textplain: One or several lines of text,  no box
##==============================================================================

textplain <- function (mid, height=0.1, lab="", adj=c(0.5, 0.5), ...)  {

  if (length (lab) == 1)
    text(mid[1],  mid[2], lab, adj=adj, ...)
  else  {
    y1      <- mid[2]+height
    ddy     <- 2*height/(length(lab)+1)
    for (i in 1:length(lab))
      text(mid[1],  y1-ddy*i, lab[i], adj=adj, ...)
  }
}
