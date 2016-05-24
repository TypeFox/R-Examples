# il fallait fields pour migraine.colors: cf remplacement splint dans spaMM.colors

`spaMM.colors` <- function (n = 64, redshift = 1) {
  orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", 
            "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", 
            "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", 
            "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", 
            "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", 
            "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
            "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", 
            "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", 
            "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
            "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", 
            "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", 
            "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
            "#AF0000", "#9F0000", "#8F0000", "#800000")
  orig[1:20] <- topo.colors(64)[1:20]
  if (n == 64 && redshift == 1) 
    return(orig)
  rgb.tim <- t(col2rgb(orig))
  temp <- matrix(NA, ncol = 3, nrow = n)
  x <- (seq(0, 1, , 64)^(redshift))
  xg <- seq(0, 1, , n)
  for (k in 1:3) {
    # hold <- splint(x, rgb.tim[, k], xg) ## using fields
    ## can be replaced by:
    colorpts <- data.frame(x=x,y=rgb.tim[,k])
    blob <- corrHLfit(y~ Matern(1|x),data=colorpts,ranFix=list(phi=1e-07,nu=4,rho=10)) ## use rho large...
    hold <- predict(blob,newdata=data.frame(x=xg))[,1] ## binding FALSE by default -> returns vector 
    ##
    hold[hold < 0] <- 0
    hold[hold > 255] <- 255
    temp[, k] <- round(hold)
  }
  rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}

## returns in canonical scale
niceLabels <- function (x, log = FALSE, lpos, maxticks = Inf, axis, base = 10L,
            ...)
  {
    if (log) {
      if (missing(lpos)) {
        lposwasmissing <- TRUE
        lpos <- c(1, 2, 5)
      }
      else lposwasmissing <- FALSE
      x <- x[x > 0]
      fc <- floor(log(min(x), base)):ceiling(log(max(x), base))
      tick <- as.vector(outer(lpos, base^fc, "*"))
      mintick <- min(tick)
      maxtick <- max(tick)
      ft <- max(c(mintick, tick[tick < min(x)]))
      lt <- min(c(maxtick, tick[tick > max(x)]))
      tick <- tick[tick >= ft/1.00000001 & tick <= lt * 1.00000001]
      if (lposwasmissing && length(tick) < 4) {
        if (axis == 1) {
          llpos <- c(1, 1.2, 1.5, 2, 3, 4, 5, 6, 8)
        }
        else llpos <- c(1, 1.2, 1.5, 2, 3, 4, 5, 6, 7, 8,
                        9)
        tick <- niceLabels(x, log = TRUE, axis = axis, lpos = llpos)
      }
      if (lposwasmissing && length(tick) > maxticks)
        tick <- niceLabels(x, log = TRUE, axis = axis, lpos = c(1,
                                                                3))
      if (lposwasmissing && length(tick) > maxticks)
        tick <- niceLabels(x, log = TRUE, axis = axis, lpos = c(1))
      if (lposwasmissing && length(tick) > maxticks) {
        base <- base * base
        tick <- niceLabels(x, log = TRUE, axis = axis, base = base,
                           lpos = c(1), maxticks = maxticks)
      }
    }
    else {
      tick <- pretty(x, high.u.bias = 0, ...)
      if (axis %in% c(1, 3) && length(tick) > 7) {
        check <- nchar(tick)
        blob <- diff(check)
        if (any(blob == 0 & check[-1] == 5)) {
          tick <- pretty(x, high.u.bias = 1.5, ...)
        }
      }
    }
    return(tick)
  }


makeTicks <- function(x, ## in canonical scale... 
                      axis, maxticks,scalefn=NULL,logticks,
                      validRange=NULL
                      ) {
  labels <- niceLabels(x, log=logticks, axis=axis) ## in canonical scale
  if( ! is.null(validRange)) {
    labels <- pmax(validRange[1],labels)
    labels <- pmin(validRange[2],labels)
  }
  if(logticks) {
    phantomat <- niceLabels(c(min(labels)/10, max(labels)*10), log=TRUE, lpos=c(1:9)) ## smalltick marks without labels
  } else phantomat <- NULL
  if( ! is.null(validRange)) {
    phantomat <- pmax(validRange[1],phantomat)
    phantomat <- pmin(validRange[2],phantomat)
  }
  at <- labels
  if( ! is.null(scalefn)) { ##not necess "log"
    at <- sapply(at,scalefn)
    invalidat <- (is.infinite(at) | is.nan(at))
    at <- at[ ! invalidat]
    labels <- labels[ ! invalidat]
    phantomat <- sapply(phantomat,scalefn)
    phantomat <- phantomat[ ! (is.infinite(phantomat) | is.nan(phantomat))]
  }
  return(list(labels=labels, at=at, phantomat=phantomat))
}

  


`spaMM.filled.contour` <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                       length.out = ncol(z)), z, xrange = range(x, finite = TRUE), 
          yrange = range(y, finite = TRUE), zrange = range(z, finite = TRUE), 
          margin=1/20,
          levels = pretty(zrange, nlevels), nlevels = 20, color.palette = spaMM.colors, 
          col = color.palette(length(levels) - 1), plot.title, plot.axes, 
          key.title=NULL, key.axes=NULL, map.asp = NULL, xaxs = "i", yaxs = "i", las = 1, 
          axes = TRUE, frame.plot = axes, ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")

  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  wmaphmap <- calc.plot.dims(x,y,xrange=xrange,yrange=yrange,margin=margin,map.asp=map.asp)  
  layout(matrix(c(2, 1), ncol = 2L), widths = c(lcm(wmaphmap[1]),lcm(wmaphmap[3])),heights=c(lcm(wmaphmap[2])),respect=TRUE)
  
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  ## SCALE
  plot.new()
  plotScale(z,levels,key.axes,key.title,axes,col)
  #
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xrange, yrange, "", xaxs = xaxs, yaxs = yaxs)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...) ## not that this function does not know 'coordinates'
  else plot.title
  invisible()
} ## spaMM.filled.contour


calc.plot.dims <- function(x,y,xrange=NULL,yrange=NULL,margin=1/20,map.asp=NULL) {
  if (is.null(xrange)) {
    xrange <- range(x)
  }
  xspan <- (xrange[2]-xrange[1])
  margex <- xspan * margin
  xrange  <- xrange+margex*c(-1,1)
  if (is.null(yrange)) {
    yrange <- range(y)
  }  
  yspan <- (yrange[2]-yrange[1])
  margey <- yspan * margin
  yrange  <- yrange+margey*c(-1,1)

  # consequences of par()
  wscale <- (3 + par("mar")[2]) * par("csi") * 2.54
  wmap <- par("din")[1]*2.54 - wscale
  Wmargin <- (par("din")[1]-par("pin")[1])*2.54
  wplotmap <- wmap - Wmargin  ## likely width of plot area
  Hmargin <- (par("din")[2]-par("pin")[2])*2.54
  max_map.asp <- ( par("din")[2]*2.54 -Hmargin)/wplotmap
  #print(paste("max_map.asp=",max_map.asp))
  
  if (is.null(map.asp)) map.asp <- yspan/xspan
  map.asp <- min(map.asp,max_map.asp)
  if (map.asp>4) map.asp <- 1
  #print(paste("map.asp=",map.asp))
  
  hmap <- wplotmap*map.asp + Hmargin
  #   if (hmap>(par("din")[2]*2.54)) {
  #     #print("hmap>(par(\"din\")[2]*2.54)")
  #     hmap <- (par("din")[2]*2.54)
  #     #reduction <- (hmap-Hmargin)/wplotmap
  #     wmap <- (hmap - Hmargin)/map.asp
  #     ## this new wmap tries to keep the aspect ratio, but it may be too narrow (and if < Wmargin, may generate a 'figure margins too large' error)
  #     if (wmap < Wmargin) {
  #       message("Aspect ratio cannot respect x and y ranges. Use 'map.asp' argument to control it directly.")
  #       wmap <- 1.05 * Wmargin ## 24/12/2014
  #     }
  #   }
  return(c(wmap,hmap,wscale))
}

plotScale <-function(z,levels,key.axes=NULL,key.title=NULL,axes,col) {
  zrange <- range(z)
  plot.window(xlim = c(0, 1), ylim = range(z),xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  
  if (is.null(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
  if (!is.null(key.title)) key.title
}

spaMMplot2D <- function (x,y,z, 
                         xrange=range(x, finite = TRUE),yrange=range(y, finite = TRUE),
                         margin=1/20,add.map= FALSE,
                         nlevels = 20, color.palette = spaMM.colors, 
                         map.asp=NULL,
                         col = color.palette(length(levels) - 1), 
                         plot.title=NULL, plot.axes=NULL, decorations=NULL,
                         key.title=NULL, key.axes=NULL, xaxs = "i", yaxs = "i", las = 1, 
                         axes = TRUE, frame.plot = axes,...) {
  dotlist <- list(...)
  ## attention,ne renvoit que ce qui correspond à l'input
  par.orig <- par(c(dotlist[intersect(names(dotlist),names(par()))],c("mar", "las", "mfrow"))) 
  on.exit(par(par.orig))
  mar.orig <- par.orig$mar
  levels <- pretty(range(z), nlevels) ## moved up to here post 1.4.4 otherwise discrepancy between main plot and scale bar 
  nlevels <- length(levels)-1
  zscaled <- 1 + floor(nlevels*(0.000001+0.999998*(z-min(z))/(max(z)-min(z)))) ## makes sure its floor( ]1,nlevels+1[ ) 
  ZColor <- color.palette(n=nlevels) ## bug corrected (spaMM.colors -> color.palette) post 1.4.4  
  wmaphmap <- calc.plot.dims(x,y,xrange=xrange,yrange=yrange,margin=margin,map.asp=map.asp)  
  layout(matrix(c(2, 1), ncol = 2L), 
         widths = c(lcm(wmaphmap[1]),lcm(wmaphmap[3])),
         heights=c(lcm(wmaphmap[2])),respect=TRUE)
  ## respect=TRUE est crucial...
  #  layout(matrix(c(2, 1), ncol = 2L), widths = c(wmap,wscale),heights=c(hmap),respect=TRUE)
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  ## SCALE
  plot.new()
  plotScale(z,levels,key.axes,key.title,axes,col)
  #
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  #  plot.new()
  ## main plot
  topontop <- order(zscaled,decreasing=FALSE) 
  plot(x=x[topontop],y=y[topontop], ## so that highest value will be printed last
       xlab="",ylab="", # give control to plot.title
       axes=FALSE, ## to retain control in later call
       xlim=xrange,ylim=yrange,xaxs = xaxs, yaxs = yaxs,
       col=ZColor[zscaled[topontop]],lwd=2)
  if (is.logical(add.map)) {
    if(add.map) {
      ## require + :: WAS the way for objects from packages in Suggests:
      ## now this is requireNamespace; but then the data such as worldMapEnv are not accessible... 
      if (requireNamespace("maps",quietly=TRUE)) {
        maps::map(,xlim=xrange,ylim=yrange,add=TRUE)  
      } else message("Package 'maps' not available, 'add.map' is ignored.")
    } 
  } else eval(add.map) ## the user may have included a map() in it but it's his problem...
  if (is.null(plot.title)) {
    do.call(title,dotlist[intersect(names(dotlist),names(formals(title)))]) ## ... may contain xlab, ylab
  } else plot.title
  if (is.null(plot.axes)) {
    if (axes) {
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  } else plot.axes
  if ( ! is.null(decorations)) decorations
  if (frame.plot) box()
  invisible()
} ## end spaMMplot2D


mapMM <- function (fitobject,Ztransf=NULL,coordinates,
                   add.points,decorations=NULL,plot.title=NULL,plot.axes=NULL,envir=-3,...) {
  ## currently add.points and decorations are equivalent:
  if ( ! missing(add.points)) warning("'add.points' is obsolete, use 'decorations'")
  if (missing(coordinates)) coordinates <- colnames(attr(fitobject,"info.uniqueGeo"))
  if (length(coordinates)!=2L) {
    stop(paste("'mapMM' plots only 2D maps, while coordinates are of length ",length(coordinates),sep=""))
  }
  pred <- predict(fitobject,binding="fitted")
  x <- pred[,coordinates[1]]
  y <- pred[,coordinates[2]]
  Zvalues <- pred[,attr(pred,"fittedName")]
  if ( ! is.null(Ztransf)) {Zvalues <- do.call(Ztransf,list(Z=Zvalues))} # 12/2014
  #dotlist <- list(...) 
  #arglist <- c(list(x=x,y=y,z=Zvalues,add.points={eval(add.points,-2)}),dotlist)
  #do.call("spaMMplot2D",arglist) # *** les eval ne passent pas dans une liste *** (PS: bc of "spaMMplot2D" instead of spaMMplot2D?)
  spaMMplot2D(x=x,y=y,z=Zvalues,
              decorations=eval(decorations,envir),
              plot.title=eval(plot.title,envir),
              plot.axes=eval(plot.axes,envir),  ## -3 -> envir in which mapMM was called.   
              ...)
}
## but list(...) has add.points=NULL for add.points originally {points.....}




`filled.mapMM` <- function(fitobject, Ztransf=NULL, coordinates, xrange = NULL, yrange = NULL, 
                           margin = 1/20, map.formula, phi = 1e-05, gridSteps = 41, 
                           decorations = quote(points(pred[, coordinates], cex = 1, lwd = 2)), 
                           add.map = FALSE, axes = TRUE, plot.axes, map.asp = NULL,
                           variance=NULL,
                           var.contour.args=list(),
                           ...) 
{
  if (missing(coordinates)) 
    coordinates <- colnames(attr(fitobject, "info.uniqueGeo"))
  if (length(coordinates) != 2L) {
    stop(paste("'map' plots only 2D maps, while coordinates are of length ", 
               length(coordinates), sep = ""))
  }
  if (length(variance)==0L) {
    variance <- list()
  } else {
    if (length(variance)>1L) stop("'variance' argument should include a single name")
    if (is.character(variance)) variance <- structure(list(TRUE),names=variance)
  }
  pred <- predict(fitobject,variances=variance, binding="fitted")
  if (missing(map.formula)) 
    map.formula <- as.formula(paste(attr(pred,"fittedName"), " ~ 1 + ", paste(findSpatial(fitobject$predictor))))
  smooobject <- corrHLfit(map.formula, data = pred, ranFix = list(phi = phi))
  smoo <- predict(smooobject,binding="dummy") ## response column does not appear to be used... 
  x <- smoo[, coordinates[1]]
  y <- smoo[, coordinates[2]]
  if (is.null(xrange)) {
    xrange <- range(x)
    margex <- (xrange[2] - xrange[1]) * margin
    xrange <- xrange + margex * c(-1, 1)
  }
  if (is.null(yrange)) {
    yrange <- range(y)
    margey <- (yrange[2] - yrange[1]) * margin
    yrange <- yrange + margey * c(-1, 1)
  }
  if (missing(plot.axes)) {
    if (axes) {
      plot.axes <- quote({title(main = "", xlab = "", ylab = "")
                         Axis(x, side = 1)
                         Axis(y, side = 2)})
    } else plot.axes <- quote(NULL)
  }
  xGrid <- seq(xrange[1], xrange[2], length.out = gridSteps)
  yGrid <- seq(yrange[1], yrange[2], length.out = gridSteps)
  newdata <- expand.grid(xGrid, yGrid)
  colnames(newdata) <- coordinates
  gridpred <- predict(smooobject, newdata = newdata,variances=variance)
  if (length(variance)==1L) {
    pvar <- attr(gridpred,names(variance))  
    varz <- matrix(pvar,ncol=length(yGrid),nrow=length(xGrid))
    contourArgs <- c(list(x=xGrid,y=yGrid,z=varz,add=TRUE),var.contour.args)
    add.varcontour <- quote(do.call(contour,contourArgs))
  } else add.varcontour <- quote(NULL)
  if (is.logical(add.map)) {
    if (add.map) {
      if (requireNamespace("maps",quietly=TRUE)) {
        add.map <- quote(maps::map(, xlim = xrange, ylim = yrange, 
                                   add = TRUE))
      }
      else {
        message("Package 'maps' not available, 'add.map' is ignored.")
        add.map <- quote(NULL)
      }
    }
    else add.map <- quote(NULL)
  }
  else add.map <- quote(add.map)
  Zvalues <- matrix(gridpred, ncol = gridSteps)
  if ( ! is.null(Ztransf)) {
    Zvalues <- do.call(Ztransf,list(Z=Zvalues))
    if (any(is.nan(Zvalues))) stop("NaN in Ztransf'ormed values: see Details of 'filled.mapMM' documentation.")
    if (any(is.infinite(Zvalues))) stop("+/-Inf in Ztransf'ormed values.")
  } 
  spaMM.filled.contour(x = xGrid, y = yGrid, z = Zvalues, margin=margin, plot.axes = {
    eval(plot.axes) ## eval in the parent envir= that of filled.mapMM; 
    ## allows ref to internal var of filled.mapMM in the call of filled.mapMM... see ?filled.mapMM 
    ## plot.new() will be evaluated before these promises are evaluated
    ## does not work if plot.new() is called one level further in a call stack...
    eval(add.varcontour)
    eval(decorations) ## evalbc it uses a local variable 'pred'
    eval(add.map)
  }, map.asp = map.asp, ...)
}
