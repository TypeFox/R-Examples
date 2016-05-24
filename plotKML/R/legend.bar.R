# Purpose        : Plot a legend bar for numeric or factor variable;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Pierre Roudier (pierre.roudier@landcare.nz); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : pre-alpha
# Note           : 1D legend with min max and avg values;
 

kml_legend.bar <- function(x, width, height, pointsize = 14, legend.file, legend.pal, z.lim = range(x, na.rm = TRUE, finite = TRUE), factor.labels, png.type = "cairo-png"){

  ## Factor-type variables:
  if(class(x) == "factor" | class(x) == "character") {
 
    z.lim <- NA
    if(missing(factor.labels)){ 
      col.no <- length(levels(as.factor(x)))  
    } else { 
      col.no <- length(factor.labels) 
    }
 
    if(missing(factor.labels)) {
      ### NOTE : This is a not a perfect implementation for a factor with a lot of categories!
      if(missing(width)) { 
        width <- max(nchar(levels(as.factor(x))))*5+70 
      } # 5 pix per character
    
      if(missing(height)) { 
        height <- length(levels(as.factor(x)))*40 
      } # 20 pix per class
    }
 
    else {
      if(missing(width)) { 
        width <- max(nchar(factor.labels))*5+70 
      } # 10 pix per character
      if(missing(height)) { 
        height <- length(factor.labels)*40 
      } # 20 pix per class
    }
 
    png(filename=legend.file, width=width, height=height, bg="transparent", pointsize=pointsize, type=png.type)
    # c(bottom, left, top, right)
    par(mar=c(.5,0,.5,1))
    plot(x=rep(1, col.no), y=1:col.no, axes=FALSE, xlab='', ylab='', pch=15, cex=4, col=legend.pal, xlim=c(0,.6*width))
  
    if(missing(factor.labels)) {
      text(x=rep(1, col.no), y=1:col.no, labels=levels(as.factor(x)), cex=.8, pos=4, offset=1, col=rgb(0.99,0.99,0.99))
    } else { 
      text(x=rep(1, col.no), y=1:col.no, labels=factor.labels, cex=.8, pos=4, offset=1, col=rgb(0.99,0.99,0.99))
    }
  
   dev.off()
  }

  ### Numeric-type variables:
  else {
    if(is.numeric(x)) {
  
    if(missing(width)) { width <- 120 } 
    if(missing(height)) { height <- 240 } 
    png(filename=legend.file, width=width, height=height, bg="transparent", pointsize=pointsize)
    par(mar=c(.5,0,.5,4))
    plot(x=0:5, y=0:5, asp=3, type="n", axes=FALSE, xlab='', ylab='')
    # get the 2-4 significant digits
    col.labels <- signif(c(z.lim[1], mean(z.lim), z.lim[2]), 2)
    color.legend(xl=0, yb=0, xr=5, yt=5, legend=col.labels, rect.col=legend.pal, gradient="y", align="rb", cex=1.4, col=rgb(0.99,0.99,0.99))
  
    dev.off()
  } 
    else { 
      stop("Vector of type 'numeric' or 'factor' expected") 
    }
  }
}

# end of script;
