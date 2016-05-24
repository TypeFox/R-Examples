##' Add color-bar to plot
##'
##' @title Add color-bar to plot
##' @param clut Color look-up table
##' @param x.range x range
##' @param y.range y range
##' @param values label values
##' @param digits number of digits
##' @param label.offset label offset
##' @param srt rotation of labels
##' @param cex text size
##' @param border border of color bar rectangles
##' @param alpha Alpha (transparency) level 0-1
##' @param position Label position left/bottom (1) or top/right (2) or no text (0)
##' @param direction horizontal or vertical color bars
##' @param \dots additional low level arguments (i.e. parsed to \code{text})
##' @export
##' @examples
##' \dontrun{
##' plotNeuro(x,roi=R,mm=-18,range=5)
##' colorbar(clut=Col(rev(rainbow(11,start=0,end=0.69)),0.5),
##'          x=c(-40,40),y.range=c(84,90),values=c(-5:5))
##'
##' colorbar(clut=Col(rev(rainbow(11,start=0,end=0.69)),0.5),
##'          x=c(-10,10),y.range=c(-100,50),values=c(-5:5),
##'          direction="vertical",border=1)
##' }
colorbar <- function(clut=Col(rev(rainbow(11,start=0,end=0.69)),alpha),
                     x.range=c(-.5,.5),y.range=c(-.1,.1),
                     values=seq(clut),digits=2,label.offset,srt=45,
                     cex=0.5,border=NA,
                     alpha=0.5,
                     position=1,
                     direction=c("horizontal","vertical"),...) {
  nlut <- length(clut)
  X <- length(agrep(tolower(direction[1]),"horizontal"))>0
  scale <- ifelse(X,diff(x.range),diff(y.range))/nlut
  barsize <- ifelse(X,diff(y.range),diff(x.range))
  if (missing(label.offset)) label.offset <- barsize/3
  delta <- ifelse(X,x.range[1],y.range[1])
  if (!is.null(values)) dM <- diff(range(values))/(nlut-1)
  for (i in seq_len(nlut+1)-1) {
      pos <- delta + (i-1)*scale
      if (X) {
          x1 <- pos; x2 <- pos+scale; y1 <- y.range[1]; y2 <- y.range[2]
      } else {
          y1 <- pos; y2 <- pos+scale; x1 <- x.range[1]; x2 <- x.range[2]
      }
      if (i>0)
          rect(x1,y1,x2,y2, col=clut[i], border=border, xpd=TRUE)
  }
  if (!is.null(values)) {
      for (i in seq_len(nlut+1)-1) {
          pos <- delta + (i-1)*scale
          rund <- format(round(min(values)+dM*i,max(1,digits)),digits=digits)
          ##      rund <- round((min(values)+dM*i)*10^digits)/(10^digits)
          x0 <- pos+(1+0.5)*scale; y0 <- y.range[2]+label.offset
          if (!X) {
              y0 <- x0;
              if (position==1) x0 <- x.range[1]-label.offset
              if (position==2) x0 <- x.range[1]+label.offset*5
              if (position==3) x0 <- x.range[1]+label.offset*1
          }
          if (i<nlut)
              text(x0,y0,rund,cex=cex,srt=srt,xpd=TRUE,...)
      }
  }
}
