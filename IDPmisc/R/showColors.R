### showColors.R

showColors <- function(col = IDPcolorRamp(20),
                       ntm = min(length(col),20),
                       border = TRUE, mar = rep(0,4))
    ## Shows Colors produced by a color vector and labels them by the
    ## index of the color vector

    ## Author: Rene Locher
    ## Version: 2009-02-19

{
  opar <- par(no.readonly = TRUE)
  par(mar = mar)
  on.exit(opar)
  lim <- 0:length(col)
  col.range <- c(0,length(col))
  plot(x=c(0,1), y=col.range,
       xlim=c(0, 1), ylim=col.range,
       type="n", ann=FALSE, axes=FALSE)

  cx <- par("cxy")[1]
  box.width <- min(3*cx,1-4*cx)
  if(border) rect(0, lim[-length(lim)], box.width,
                  lim[-1], col=col) else
  rect(0, lim[-length(lim)], box.width, lim[-1], col=col, border=col)
  if(ntm>0){
    ap <- pretty(lim[-1],n=ntm)
    text(box.width+cx, ap-0.5, labels=paste(ap),adj=0)
  }
} ## showColors
