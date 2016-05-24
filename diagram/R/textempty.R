
##==============================================================================
# textempty: One or several lines of text,  background colored,  no box
##==============================================================================

textempty <- function (mid, lab="", adj=c(0.5, 0.5),
    box.col="white", cex=1, ...)   {

  text.width  <- max(strwidth(lab,  units = "user",  cex = cex))
  text.height <- strheight(lab,  units = "user",  cex = cex)

  rect(mid[1]-text.width*adj[1],  mid[2]-text.height*adj[2],
       mid[1]+text.width*(1-adj[1]),  mid[2]+text.height*(1-adj[2]),
       col=box.col, border=NA)

  textplain(mid=mid, height=text.height, lab=lab, adj=adj, cex=cex, ...)

}
