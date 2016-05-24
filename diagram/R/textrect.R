
##==============================================================================
# textrect: One or several lines of text,  in a rectangle
##==============================================================================

textrect <- function (mid, radx, rady=radx*length(lab), lwd=1,
    shadow.size=0.01, adj=c(0.5, 0.5), lab="", box.col="white",
    lcol="black", shadow.col="grey", angle=0, ...)  {

  shadowbox("rect", mid=mid, radx=radx, rady=rady,
            shadow.size=shadow.size, shadow.col=shadow.col, box.col=box.col,
            lcol=lcol, lwd=lwd, angle=angle)
  textplain(mid, rady, lab, adj, ...)

}

textparallel <- function (mid, radx, rady=radx*length(lab), lwd=1,
    shadow.size=0.01, adj=c(0.5, 0.5), lab="", box.col="white",
    lcol="black", shadow.col="grey", angle=0, theta = 90, ...)  {

  shadowbox("parallel", mid=mid, radx=radx, rady=rady,
            shadow.size=shadow.size, shadow.col=shadow.col, box.col=box.col,
            lcol=lcol, lwd=lwd, angle=angle, theta = theta)
  textplain(mid, rady, lab, adj, ...)
}
