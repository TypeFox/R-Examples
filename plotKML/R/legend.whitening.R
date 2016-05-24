# Purpose        : Produce a PNG legend file for whitening;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : pre-alpha
# Note           : this technique requires a special 2D legend;

kml_legend.whitening <- function(legend.res = 0.01, width=120, height=300, pointsize = 14, x.lim, e.lim, leg.asp = 0.3*width/height, legend.file = "whitening_legend.png", matte = FALSE, png.type="cairo-png"){
  
  xlg <- seq(.01,1,by=legend.res)
  ylg <- seq(.01,1,by=legend.res)
  # empty grid
  leg <- expand.grid(xlg, ylg, KEEP.OUT.ATTRS=FALSE)
  # Hues
  f1 <- -90-leg[,2]*300
  f2 <- ifelse(f1<=-360, f1+360, f1)
  H <- ifelse(f2>=0, f2, (f2+360))
  # Saturation
  S <- 1-leg[,1]
  # Intensity
  V <- 0.5+leg[,1]/2
  HSV <- as.vector(t(matrix(hex(HSV(rev(H), S, V)), nrow=length(ylg), ncol=length(xlg))))
  leg.plt <- pixmapIndexed(data=1:length(HSV), nrow=length(ylg), ncol=length(xlg), bbox=c(e.lim[1], x.lim[1], e.lim[2], x.lim[2]), col=HSV)
  # par(las = 0)
  
  png(filename=legend.file, width=width, height=height, bg="transparent", pointsize=pointsize, type=png.type)
  par(mar=c(2.5,2.5,0.5,0))
  plot(leg.plt, axes=FALSE, col.lab=rgb(0.99,0.99,0.99), bg=NA, asp=leg.asp)
  axis(side=1, at=e.lim, cex=.8, col.axis=rgb(0.99,0.99,0.99), col.lab=rgb(0.99,0.99,0.99))
  axis(side=2, at=signif(x.lim, 3), cex=.8, col.axis=rgb(0.99,0.99,0.99), col.lab=rgb(0.99,0.99,0.99))
  dev.off()
 
  ## Force transparency (requires ImageMagick):
  if(matte==TRUE){
  convert <- get("convert", envir = plotKML.opts)
  if(nchar(convert)==0){
    plotKML.env(silent = FALSE, show.env = FALSE)
    convert <- get("convert", envir = plotKML.opts)
    if(!nchar(convert)==0){
      system(paste(convert, ' ', legend.file, ' -matte -transparent "#FFFFFF" ', legend.file, sep=""))
    }
  }
  else { 
    warning("PNG transparency possibly ineffective. Install ImageMagick and add to PATH. See '?kml_legend.whitening' for more info.")
  } 
  }  
   
}

# end of script;
