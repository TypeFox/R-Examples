
##==============================================================================
# texthexa: One or several lines of text,  in a hexa thingie
##==============================================================================

texthexa <- function (mid, radx, rady=radx*length(lab), lwd=1,
    shadow.size=0.01, adj=c(0.5, 0.5), lab="", box.col="white",
    lcol="black", shadow.col="grey", angle=0, ...)  {

  shadowbox("hexa", mid=mid, radx=radx, rady=rady,
            shadow.size=shadow.size, shadow.col=shadow.col, box.col=box.col,
            lcol=lcol, lwd=lwd, angle=angle)
  textplain(mid=mid, height=rady, lab=lab, adj=adj, ...)

}
