
##==============================================================================
# textround: One or several lines of text,  in a rectangle with rounded sides
##==============================================================================

textround <- function (mid, radx, rady=radx*length(lab), lwd=1,
    shadow.size=0.01, adj=c(0.5, 0.5), lab="", box.col="white",
    lcol="black", shadow.col="grey", angle=0, rx=rady, ...)  {

  shadowbox("round", mid=mid, radx=radx, rady=rady,
            shadow.size=shadow.size, shadow.col=shadow.col, box.col=box.col,
            lcol=lcol, lwd=lwd, rx=rx, angle=angle)
  textplain(mid, rady, lab, adj, ...)

}
