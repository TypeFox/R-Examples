#' Write text to plot with halo underneath
#' 
#' Write text to plot. A field the size of each label is drawn beneath it, so
#' the text can be read easily even if there are many points in the plot.
#' Fields can be rectangular, elliptic or rectangular with roundeed edges.
#' 
#' @details Specifying pos and offset will currently change the position of the text, but not of the field.\cr 
#' srt is not supported yet.\cr 
#' lend, ljoin and lmitre can not be specified for rect, to keep argument number low.\cr 
#' density (crosshatch etc.) is not supported, as this would distract from the text.
#' # Search Engine Keywords:
#' R Text visible on top
#' R labeling with color underneath
#' R Creating text with a halo
#' R Text with shadow
#' 
#' @return None
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, April 2013 + March 2014
#' @seealso \code{\link{text}}; \code{shadowtext} in package \code{TeachingDemos}, see \url{http://stackoverflow.com/questions/25631216};
#'          \code{s.label} in package \code{ade4}, which is not so versatile and doesn't work with logarithmic axes
#' @references with inspiration taken from \code{ordilabel} in package \code{vegan} and thanks to Jari Oksanen for his comments
#' @keywords aplot
#' @export
#' @examples
#' 
#' # TextFields with mixed field shapes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' set.seed(13);  plot(cumsum(rnorm(100)), type="l", main="berryFunctions::textField")
#' for(i in 2:7) lines(cumsum(rnorm(100)), col=i)
#' textField(40, 4, "default")
#' textField(40, 0, "some options", col=2, fill=4, margin=c(-0.4, 0.9), font=2)
#' # Ellipsis (looks better in vector graphics like pdf):
#' textField(80, 2, "field='ellipse'", field="ell", mar=c(0.5, 2.3), border=5)
#' # Rectangular field with edges rounded:
#' textField(60,-3, "field='Rounded'", field="rounded", fill="orange", cex=1.7)
#' 
#' # Field type can be abbreviated (partial matching), margin may need adjustment:
#' textField(90, 5, "short", field="ell", fill=7, border=4, mar=-0.4)
#' 
#' # Rounded can also vectorized:
#' textField(30, c(2,0,-2,-4,-6), paste("rounding =",seq(0,1,len=5)), field="round",
#'     fill=(2:6), mar=1, rounding=seq(0,1,len=5), border=1)
#' # turn off warning about recycling:
#' textField(80, c(-5,-6.5), c("Ja", "Nein"), field="round", fill=6:8, quiet=TRUE)
#' 
#' 
#' set.seed(007); plot(rnorm(1e4)) ; abline(v=0:5*2e3, col=8)
#' # Default settings:
#' textField(5000, 0, "Here's some good text")
#' # right-adjusted text (the field box still extends 'margin' stringwidths em):
#' textField(2000, -1, "Some more (smores!)", cex=1.5, adj=0, col=2)
#' # Field color, no extra margin beyond baseline (excluding descenders):
#' textField(2000, -2, "more yet", col=2, fill="blue", margin=0)
#' # margin can be one number for both x and y direction ... :
#' textField(1000, 2, "Up we go", fill=7, margin=1.4)
#' # ... or two (x and y different), even negative:
#' textField(5000, 2, "to the right", col=2, fill=4, margin=c(-0.4, 0.9))
#' # Fonts can be set as well:
#' textField(5000, 1, "And boldly down in bold font", font=2, border=3)
#' # Text can expand outsinde of the plot region (figure) into the margins:
#' textField(11000, -2, "Hi, I'm a long block of text", adj=1, fill="red")
#' textField(11000, -3, "You're not outside the plot!", adj=1, xpd=TRUE, fill="red")
#' # And most parameters can be vectorized, while x/y are recycled:
#' textField(3000, c(-3, -3.7), c("0", "good"), border=c("red",3), lty=1:2)
#' 
#' 
#' # textField even works on logarithmic axes:
#' mylabel <- c("This","is (g)","the","ever-\n great","Sparta")
#' plot(10^runif(5000, -1,2), log="y", col=8)
#' textField(1000, c(100,20,4,2,0.5), mylabel, fill=2, mar=0, expression=FALSE)
#' textField(2500, c(100,20,4,2,0.5), mylabel, fill=4, mar=0, expression=TRUE)
#' textField(4000, c(100,20,4,2,0.5), mylabel, fill=3, mar=0)
#' textField(c(1,2.5,4)*1000, 0.2, paste("expression=\n", c("FALSE","TRUE","NA")))
#' 
#' # In most devices, vertical adjustment is slightly off when the character string
#' # contains no descenders. The default is for centered text:  adj = c(0.5, NA).
#' # For drawing the field, adj[2] is in this case set to 0.5.
#' # Text positioning is different for NA than for 0.5, see details of ?text
#' # I'm working on it through expression, which does not work with newlines yet
#' 
#' @param x X coordinates, if necessary, they are recycled
#' @param y Y coordinates
#' @param labels labels to be placed at the coordinates, as in \code{\link{text}}. DEFAULT: seq_along(x)
#' @param fill fill is recycled if necessary. With a message when quiet = FALSE. DEFAULT: "white"
#' @param border ditto for border. DEFAULT: NA
#' @param expression If TRUE, labels are converted to expression for better field positioning through expression bounding boxes.  
#'        If NA, it is set to TRUE for labels without line breaks (Newlines, "\\n").
#'        If FALSE, no conversion happens. DEFAULT: NA
#' @param margin added field space around words (multiple of em/ex). DEFAULT: 0.3
#' @param field 'rectangle', 'ellipse', or 'rounded', partial matching is performed. DEFAULT: "rect"
#' @param nv number of vertices for field = "ellipse" or "rounded". low: fast drawing.
#'        high: high resolution in vector graphics as pdf possible. DEFAULT: 1000
#' @param rounding between 0 and 1: portion of height that is cut off rounded at edges when field = "rounded". DEFAULT: 0.75
#' @param lty line type. DEFAULT: par("lty")
#' @param lwd line width. DEFAULT: par("lwd")
#' @param cex character expansion. DEFAULT: par("cex")
#' @param xpd expand text outside of plot region ("figure")?. DEFAULT: par("xpd")
#' @param adj vector of length one or two. DEFAULT: par("adj")
#' @param pos in 'text', pos overrides adj values. DEFAULT: NULL
#' @param offset I want the field to still be drawn with adj, but have it based on pos. DEFAULT: 0.5
#' @param quiet Suppress warning when Arguments are recycled?. DEFAULT: FALSE
#' @param \dots further arguments passed to strwidth and text, like font, vfont, family
#' 
textField <- function(
x,
y,
labels=seq_along(x),
fill="white",
border=NA,
expression=NA,
margin=0.3,
field="rect",
nv=1000,
rounding=0.75,
lty=par("lty"),
lwd=par("lwd"),
cex=par("cex"),
xpd=par("xpd"),
adj=par("adj"),
pos=NULL,
offset=0.5,
quiet=FALSE,
...)
{
# Partial matching field--------------------------------------------------------
PossibleValues <- c("rectangle", "ellipse", "rounded")
field <- PossibleValues[pmatch(field,  PossibleValues)]
#
# Recycling --------------------------------------------------------------------
# Recycle x or y, if one is shorter than the other. Code taken from xy.coords
nx <- length(x)  ;  ny <- length(y)
if( nx < ny )
  { x <- rep(x, length.out=ny)  ; nx <- length(x) }
  else
  y <- rep(y, length.out=nx)
if(length(labels) > nx) stop("more labels than coordinates were given.")
# Recycle arguments if necessary
if(! quiet)
 {
 rarg <- NULL # recycled arguments
 rtim <- NULL # number of times recycled
 if(length(labels)!=nx){rarg<-c(rarg,"labels"); rtim<-c(rtim,nx/length(labels))}
 if(length(fill)  !=nx){rarg<-c(rarg,"fill")  ; rtim<-c(rtim,nx/length(fill))}
 if(length(border)!=nx){rarg<-c(rarg,"border"); rtim<-c(rtim,nx/length(border))}
 if(length(lty)   !=nx){rarg<-c(rarg,"lty")   ; rtim<-c(rtim,nx/length(lty))}
 if(length(lwd)   !=nx){rarg<-c(rarg,"lwd")   ; rtim<-c(rtim,nx/length(lwd))}
 if(length(cex)   !=nx){rarg<-c(rarg,"cex")   ; rtim<-c(rtim,nx/length(cex))}
 if(length(rounding)!=nx & field=="rounded"){rarg<-c(rarg,"rounding");rtim<-c(rtim,nx/length(rounding))}
 if(!is.null(rarg))
  {
  rtim <- round(rtim, 2)
  warning("The following arguments have been recycled:\n",
          paste(format(rarg), format(rtim), "times", collapse="\n"))
  }                         # formatC(rtim, width=5)
 }# end if !quiet
if(length(labels) != nx ) labels <- rep(labels, length.out=nx)
if(length(fill)   != nx )   fill <- rep(  fill, length.out=nx)
if(length(border) != nx ) border <- rep(border, length.out=nx)
if(length(lty)    != nx )    lty <- rep(   lty, length.out=nx)
if(length(lwd)    != nx )    lwd <- rep(   lwd, length.out=nx)
if(length(cex)    != nx )    cex <- rep(   cex, length.out=nx)
if(length(rounding) != nx & field=="rounded") rounding <- rep(rounding, length.out=nx)
#
# Dimensioning -----------------------------------------------------------------
# better field positioning through expression bounding boxes:
 labels2 <- as.list(labels)
# labels without newline:
nl <- which(!sapply(labels2, grepl, pattern="\n"))
if(is.na(expression) & length(nl)>0 )
  for(i in nl) labels2[[i]] <- as.expression(labels2[[i]])
if(isTRUE(expression)) labels2 <- lapply(labels2, as.expression)
#
# Dimension of the character string in plot units:
w <- sapply(labels2, strwidth , cex=cex, ...)
h <- sapply(labels2, strheight, cex=cex, ...)
# Box height times number of line breaks
if(isTRUE(expression))
##   h <- h + h*sapply(gregexpr("\n", labels), function(x) sum(x>0)) # false
##labels=c("Bug","oo-\nbahg", "Bug-\nbahg\ngrh")
h <- sapply(as.list(labels), function(xx) 
    {
    xxsplit <- strsplit(xx, "\n")[[1]]
    sum(sapply( lapply(as.list(xxsplit), as.expression), strheight, cex=cex, ...))
    })
#browser() # this sometimes is a list!



#h <- strheight("lg", cex=cex, ...)
# Extra-space (margins) around characters:
if(field=="ellipse") margin <- margin + 1.5 # bigger margin for ellipses needed
if(is.na(margin[2])) margin[2] <- margin[1]
mar_x <- strwidth ("m", cex=cex, ...)*margin[1]  # 'em' = stringwidth of letter m
mar_y <- strheight("x", cex=cex, ...)*margin[2]  # in typography, analog for ex
# For logarithmic axes:
if(par("ylog")) y <- log10(y)
if(par("xlog")) x <- log10(x)
# Adjust the adj parameter for vertical positioning (if only one value is given):
if(is.na(adj[2]))
  adjy <- 0.5
  else
  adjy <- adj[2]
# Adjust adj parameter based on pos and offset:
# if(!is.null(pos))
#
# Drawing ----------------------------------------------------------------------
if(field=="rectangle")
# Plot rectangular fields: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{ xleft <- x-w*adj[1]-mar_x;  xright <- x-w*adj[1]+w+mar_x
ybottom <- y-h*adjy  -mar_y;    ytop <- y-h*adjy  +h+mar_y
rect(  xleft = if(par("xlog")) 10^xleft   else xleft,
      xright = if(par("xlog")) 10^xright  else xright,
     ybottom = if(par("ylog")) 10^ybottom else ybottom,
        ytop = if(par("ylog")) 10^ytop    else ytop,
        col=fill, border=border, xpd=xpd, lty=lty, lwd=lwd)
}else if(field=="ellipse")
#
# Plot elliptic fields: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
for(i in 1:nx)
  {xell <- x[i] + (w[i]+mar_x)/2*cos(seq(0, 2*pi, length=nv))
   yell <- y[i] + (h[i]+mar_y)/2*sin(seq(0, 2*pi, length=nv))
   polygon(x = if(par("xlog")) 10^xell else xell,
           y = if(par("ylog")) 10^yell else yell,
           col=fill[i], border=border[i], xpd=xpd, lty=lty[i], lwd=lwd[i])
  }
}
else if(field=="rounded")
#
# Plot rectangular fields with rounded corners: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
#if(par("ylog") | par("xlog")) stop("Rounded not yet possible with log axes")
if(any(rounding < 0 || rounding > 1  )) stop("Wrong rounding value. Needs to be between 0 and 1.")
asp <- diff(par("usr")[3:4])/diff(par("usr")[1:2]) # current aspect ratio y to x
for(i in 1:nx)
  {
  XL <- x[i] - w[i]*adj[1] - mar_x        # Left x position, as with rectangle
  XR <- x[i] - w[i]*adj[1] + w[i] + mar_x # right x
  YB <- y[i] - h[i]*adjy   - mar_y        # bottom y position
  YT <- y[i] - h[i]*adjy   + h[i] + mar_y # top y
  xi <- rounding[i]*(mar_x+h[i]/asp/2)  # x inset of rounded corner
  yi <- rounding[i]*(mar_y+h[i]/2)      # y inset
  elx <- function(from,to) xi*cos(seq(from,to,length.out=nv/4)) # elliptic corners function
  ely <- function(from,to) yi*sin(seq(from,to,length.out=nv/4))
  # x and y coordinates:
  xc <- c(XR-xi+elx(0,pi/2), XR-xi, XL+xi, XL+xi+elx(pi/2,pi), XL,    XL,
   XL+xi+elx(pi,3*pi/2), XL+xi, XR-xi, XR-xi+elx(3*pi/2,2*pi), XR,    XR)
  yc <- c(YT-yi+ely(0,pi/2), YT,    YT,    YT-yi+ely(pi/2,pi), YT-yi, YB+yi,
   YB+yi+ely(pi,3*pi/2), YB,    YB,    YB+yi+ely(3*pi/2,2*pi), YB+yi, YT-yi)
  polygon(x = if(par("xlog")) 10^xc else xc,
          y = if(par("ylog")) 10^yc else yc,
          col=fill[i], border=border[i], xpd=xpd, lty=lty[i], lwd=lwd[i])
  } # End of for loop
}
else stop("Wrong field type specified. Use 'rectangle', 'ellipse', or 'rounded'.")
#
# Writing ----------------------------------------------------------------------
# Write text:
text(if(par("xlog")) 10^x else x,  if(par("ylog")) 10^y else y,
     labels=labels, cex=cex, xpd=xpd, adj=adj, pos=pos, offset=offset, ...)
} # End of function

