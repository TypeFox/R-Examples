
##==============================================================================
# textdiamond: One or several lines of text,  in a diamond
##==============================================================================

textdiamond <- function (mid, radx, rady=NULL, lwd=1, shadow.size=0.01,
    adj=c(0.5, 0.5),lab="", box.col="white", lcol="black",
    shadow.col="grey", angle=0, ...)    {

  if (is.null(rady)) {
    pin <- par("pin")
    rady  <- radx*pin[1]/pin[2]* length(lab)
  }

  shadowbox("diamond", mid=mid, radx=radx, rady=rady,
            shadow.size=shadow.size, shadow.col=shadow.col, box.col=box.col,
            lcol=lcol, lwd=lwd, angle=angle)
  textplain(mid, rady, lab, adj, ...)

}
