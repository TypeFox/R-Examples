.packageName <- "hydrogeo"
# Copyright Myles English 2008-2014
# Contact: myles@rockhead.biz
# License: BSD_2_clause

# Definition of piper class
# If cl is used instead of Cl then get
# "Error in .getClassFromCache(Class, where) :
#    Class should be either a character-string name or a class definition"
setClass("piper",
         representation(size="numeric",
                        Ca="vector",
                        Mg="vector",
                        Cl="vector",
                        SO4="vector",
                        WaterType="vector",
                        group="vector",     # TODO: make a factor
                        IDs="vector",       # NOTE: row.names changed to IDs
                        pt.col="vector",    # colours for WaterType
                        pt.pch="vector",    # symbols for WaterType
                        call="call"    # call that created it
                        ),
         prototype(size=300)
         )

##_______________INITIALISER ____________________
# @WaterType should be the first factor
# @group - should be the second factor, to allow plotting on different plots
# when there are a lot of points but the colours and symbols need to be
# preserved between plots.
# i.e. col and pch by @WaterType
#      plot by @group
# WaterTypes is shown by colour
# Sample ID is shown by colour and pch via legend with IDs

setMethod(f="initialize",
          signature="piper",
          definition=function(.Object, l, call=sys.call(), size=NULL,
              group=NULL, colours=NULL, numbersymbols=FALSE,
              wt.col=NULL, wt.pch=NULL, pt.col=NULL)
          {
            # initialise member data
            if ( !is.null(size) ) {
              .Object@size <- size
            } else {
                .Object@size <- 1
            }
            .Object@Ca <- l$Ca
            .Object@Mg <- l$Mg
            .Object@Cl <- l$Cl
            .Object@SO4 <- l$SO4
            .Object@pt.pch <- l$pt.pch
            if ( !is.null(l$group) ) {
                .Object@group <- group
            } else { # plot all samples on the same plot
                .Object@group <- rep( 1, times=length(l$Ca) )
            }
            if ( !is.null(l$WaterType) ) {
                .Object@WaterType <- l$WaterType
            } else { # treat each sample as an individual water type
                .Object@WaterType <- seq( 1:length(row.names(l)) )
            }
            if ( !is.null(l$IDs) ) {
                .Object@IDs <- l$IDs
            } else {
                .Object@IDs <- row.names(l)
            }
            # TODO: make this less confusing
            # make the colours and symbols for the points now
            # because it saves confusion later
            # wt.col and wt.pch should be factors (?)
            # pt.col and pt.pch should be vectors
            # pt. overides wt.
          
            # set pt.col
            if ( ! is.null( l$pt.col ) ) {              # if specified
              .Object@pt.col <- l$pt.col                # assign it
            } else {                                    # else calculate it...
              wtf <- as.factor(.Object@WaterType)
              if ( ! is.null( l$wt.col ) ) {
                if ( length(l$wt.col) != length( levels(wtf) ) ) {
                  cat("ERROR: wt.col wrong length for WaterType!")
                  return(invisible())
                } else { levels(wtf) <- l$wt.col }
              } else {  # the default for pt.col
                levels( wtf ) <- seq( 1:length(levels(wtf)) )
              }
                  
              .Object@pt.col <- as.vector( wtf )
            }

            # set pt.pch
            # input wt:      2 2 1 2 3
            # output pt.col  2 2 1 2 3
            #        pt.pch  1 2 1 3 1
            if ( ! is.null( l$pt.pch ) ) {             # if specified
              .Object@pt.pch <- l$pt.pch               # assign it
            } else {                                   # else calculate it...
              wtf <- as.factor(.Object@WaterType)
              pch <- .Object@WaterType                 # initialise
              ## if ( ! is.null( l$wt.pch ) ) {
##               for ( i in levels(wtf) ) {
##                  # get subset , replace with values from vector
##                  lrp<-sum( wtf==i ) # the number of samples of that watertype
##                  pch[ wtf==i ] <- l$wt.pch[1:lrp]
##                }
##              } else {
##                # loop through levels
##                for ( i in levels(wtf) ) {
##                  lrp<-sum( wtf==i ) # the number of samples of that watertype
##                  pch[ wtf==i ] <- seq(lrp)
##                }
##              }
           ###   .Object@pt.pch <- pch                 # assign
          }
            .Object@call <- call
              return(.Object)
          }
)

# Sets up the plot paper with two triangles and a diamond
plotpaper =
  function(x, main=NULL, xyaxes=FALSE, ...)
  {
    p <- (x@size/11)
    r <- (x@size/22)

    plot.default(0, 0, type="n", axes=xyaxes, lty=1, lwd=1,
                 xlim=c(0, x@size + p), ylim=c(-p, x@size),
                 frame.plot=FALSE, ann=TRUE, ylab="", xlab="", ...)
    if(is.null(main))
        main <- x@call
    title(main=main)

    # draws grid lines
    thickxf <- c(0, (10 * r), (5 * r), (12 * r), x@size, (17 * r),
                 (x@size/2), (16 * r), (x@size/2), (6 * r))
    thickxt <- c((10 * r), (5 * r), 0, x@size, (17 * r), (12 * r),
                 (16 * r), (x@size/2), (6 * r), (x@size/2))
    thickyf <- c(0, 0, (10 * r), 0, 0, (10 * r), (2 * r), (12 * r), 
                 x@size, (12 * r))
    thickyt <- c(0, (10 * r), 0, 0, (10 * r), 0, (12 * r), x@size, 
                 (12 * r), (2 * r) )
    
    xf <- c(thickxf, (2 * r), (4 * r), (6 * r), (8 * r), (14 * r), 
            (16 * r), (18 * r), (20 * r), (21 * r), (20 * r), (19 * r),
            (18 * r), (21 * r), (20 * r), (19 * r), (18 * r), (9 * r), (8 * r),
            (7 * r), (6 * r), (9 * r), (8 * r), (7 * r), (6 * r), (7 * r),
            (8 * r), (9 * r), (10 * r), (7 * r), (8 * r), (9 * r), (10 * r))
    
    xt <- c(thickxt, r, (2 * r), (3 * r), (4 * r), (13 * r), (14 * r), (15 * r),
            (16 * r), (13 * r), (14 * r), (15 * r), (16 * r), (20 * r),
            (18 * r), (16 * r), (14 * r), (8 * r), (6 * r), (4 * r), (2 * r),
            r, (2 * r), (3 * r), (4 * r), (12 * r), (13 * r), (14 * r),
            (15 * r), (12 * r), (13 * r), (14 * r), (15 * r))
    
    yf <- c(thickyf, 0, 0, 0, 0, 0, 0, 0, 0, (2 * r), (4 * r), (6 * r), (8 * r),
            (2 * r), (4 * r), (6 * r), (8 * r), (2 * r), (4 * r), (6 * r),
            (8 * r), (2 * r), (4 * r), (6 * r), (8 * r), (14 * r), (16 * r),
            (18 * r), (20 * r), (10 * r), (8 * r), (6 * r), (4 * r))
    
    yt <- c(thickyt, (2 * r), (4 * r), (6 * r), (8 * r), (2 * r), (4 * r),
            (6 * r), (8 * r), (2 * r), (4 * r), (6 * r), (8 * r), 0, 0, 0, 0, 0,
            0, 0, 0, (2 * r), (4 * r), (6 * r), (8 * r), (4 * r), (6 * r),
            (8 * r), (10 * r), (20 * r), (18 * r), (16 * r), (14 * r))

    segments(xf, yf, xt, yt, col=par("fg"), lty=1, lwd=par("lwd"))
  }



# Add axes to the plot TODO: make this more flexible
setMethod(
          f="Axis",
          signature=c("piper"),
          definition=function(x=NULL, at=NULL, ..., side,
                              labels=NULL)
   {
    p <- (x@size/11)      # for scaling
    r <- (x@size/22)      # for scaling

    if ( is.null(labels) ) { # hack. TODO Add comment why this is necessary
      addLabels <- TRUE
      labels <- NULL
    }
    cex.axis=0.7 # FIXME
    # add axis titles
    vfont <-  c("serif", "italic")
    
    # label bottom axes (the two horizontal axes)
    xstr <- c(5 * r, 17 * r)
    ystr <- c(-r, -r)

    text(xstr, ystr, labels=c( expression(Ca^{2+''}), expression(paste(Cl[2],
                               ''^{symbol('-')})) ),
         vfont=vfont, cex=cex.axis)

    # label axes parallel with a line angled at 60% from horizontal
    # in a clockwise direction
    xgh <- c(14 * r, 8 * r, 20 * r)
    ygh <- c(18 * r, 6 * r, 6 * r)

    text(xgh, ygh,
         labels=c( expression( paste(Ca^{2+''}, " ", symbol('&'), " ",
                                       Mg^{2+''})),
                     expression( paste(Na^{symbol('+')}, " ", symbol('&'), " ",
                                       K^{symbol('+')})),
                     expression( paste(SO[4],''^{2-''}))),
         srt=300, vfont=vfont, cex=cex.axis)

    # label axes parallel with a line angled at 60% (approx.) from horizontal
    # in an anti clockwise direction    
    xla <- c(14 * r, 8 * r, 2 * r)
    yla <- c(6 * r, 18 * r, 6 * r)
    
    text(xla, yla,
         labels=c( expression( paste(CO[3],''^{2-''}, " ", symbol('&'), " ",
                                 HCO[3],''^{2-''})),
                     expression( paste(SO[4],''^{2-''}, " ", symbol('&'), " ",
                                 Cl[2], ''^{symbol('-')})),
                     expression( paste(Mg^{2+''}))),
         srt=60, vfont=vfont, cex=cex.axis)

    if ( addLabels != FALSE ) labelAxes(x)
  }
  )

setGeneric (
    "labelAxes",
    function(x, cex.lab=0.1, side=-1, ...)
    standardGeneric("labelAxes")
    )

# was LabelAxes
setMethod(
          f="labelAxes",
          signature="piper",
          definition=function(x, cex.lab=0.35, side=-1, ...)
  {
    p <- (x@size/11)      # for scaling
    r <- (x@size/22)      # for scaling

    #TODO: Most of these don't work

    ## Tick mark labels:
    # add axis tick mark labels: 20, 40, 60, 80 percent labels
    adj <- 0.5     # to centre text on the point
    labels <- c(20 * (1:4))  # too clever by half: 20, 40, 60, 80 percent
    vfont <- c("serif", "italic")

    # for calculating offsets from line intersections
    ofs <- function(deg) {
      dd <- tan( (deg/360) * 2 * pi )
      offset <- cex.lab * dd * (x@size/50)
      return( offset )
    }

    if(side==-1 || side==1){
    # flat labels LHS - Mg2+
    xe <- c(r * 1:4)
    ye <- c(p * 1:4) 
    text(xe, ye, labels=labels, pos=2, offset=0.1, vfont=vfont, cex=cex.lab)
}
    if(side==-1 || side==2){
    # flat labels RHS - SO42-
    xe <- c(r * 21:18)
    ye <- c(p * 1:4) 
    text(xe, ye, labels=labels, pos=4, offset=0.1, vfont=vfont, cex=cex.lab)
}
    if(side==-1  || side==3){
    # labels rotated 60degrees
    # Ca2+ & Mg2+
    srt <- 60
    delta <- ofs(srt)
    xd <- c(r * 15:12) + delta
    yd <- c(p * 7:10) + delta
    text(xd, yd, labels=labels, vfont=vfont, srt=srt, cex=cex.lab, adj=adj)
}
    if(side==-1 || side==4){
    # labels rotated 120degrees
    # Ca2+
    srt <- 120
    delta <- ofs(srt)
    xa <- c(p * 4:1) - delta
    ya <- c(0) + delta
    text(xa, ya, labels=labels, vfont=vfont, srt=srt, cex=cex.lab, adj=adj)
}

    if(side==-1 || side==5){
    # CO32- & HCO32-, triangle
    srt <- 120
    delta <- ofs(srt)
    xa <- c(r * 16:13) + delta
    ya <- c(p * 4:1) - delta
    text(xa, ya, labels=labels, vfont=vfont, srt=srt, cex=cex.lab, adj=adj)
}
    
    if(side==-1 || side==6){
    # CO32- & HCO32-, diamond
    srt <- 120
    delta <- ofs(srt)
    xa <- c(r * 15:12) - delta
    ya <- c(p * 5:2) + delta
    text(xa, ya, labels=labels, vfont=vfont, srt=srt, cex=cex.lab, adj=adj)
}
    if(side==-1 || side==7){
    # labels rotated 240degrees
    # Cl-
    srt <- 240
    delta <- ofs(srt)
    xb <- c(p * 7:10) - delta
    yb <- c(0, 0, 0, 0) - delta
    text(xb, yb, labels=labels, vfont=vfont, srt=srt, cex=cex.lab, adj=adj)
}
    if(side==-1 || side==8){
    # Na+ & K+, triangle
    srt <- 240
    delta <- ofs(srt)
    xb <- c(r * 6:9) + delta
    yb <- c(p * 4:1) + delta
    text(xb, yb, labels=labels, vfont=vfont, srt=srt, cex=cex.lab, adj=adj)
}
    if(side==-1 || side==9){
    # Na+ & K+, diamond
    srt <- 240
    delta <- ofs(srt)
    xb <- c(r * 7:10) - delta
    yb <- c(p * 5:2) - delta
    text(xb, yb, labels=labels, vfont=vfont, srt=srt, cex=cex.lab, adj=adj)
}
    if(side==-1 || side==10){
    # labels rotated 300degrees
    # SO42- & Cl-
    srt <- 300
    delta <- ofs(srt)
    xc <- c(r * 7:10) + delta
    yc <- c(p * 7:10) - delta
    text(xc, yc, labels=labels, vfont=vfont, srt=srt, cex=cex.lab, adj=adj)
}
  }
  )


setMethod(
          f="points",
          signature="piper",
          definition=function(o, x, col=c(1), cex=0.75, ...)
    {
    # get coordinates for cation triangles  FIXME: length(x[0,])/1100)
    cationy <- c(x$Mg * 5 * o@size/1100)
    cationx <- c((5 * o@size/11) * (1 - (x$Ca/100) + (x$Mg/200)))

    # get coordinates for anion triangles
    aniony <- c(x$SO4 * 5 * o@size/1100)
    anionx <- c((6 * o@size/11) + ((5 * o@size/11) * (x$Cl/100)) + (1/2) *
                (5 * o@size/11) * (x$SO4/100))

    # get coordinates for points projected on the diamond
    projx <- ( (1/2) * (cationx + anionx) + ( (1/4) * (aniony - cationy) ) )
    projy <- ( anionx - cationx + ( (1/2) * (aniony + cationy) ) )

    px <- c(anionx, cationx, projx)
    py <- c(aniony, cationy, projy)

    # plot all the points in one go
    points.default(px, py, type="p", lty=1, lwd=1, pch=x$pt.pch,
                   col=x$pt.col, cex=cex, ...)
  }
  )


setMethod(
          f="plot",
          signature="piper",
          definition=function(x, y, ..., axes=TRUE, points=TRUE, cex.axis=0.7,
                              cex.lab=0.5, xyaxes=FALSE, side=-1)
                                        # group=NULL,
          {
            plotpaper(x, xyaxes=xyaxes, ...)
            # plot axes
            #print("side in plot() is:")
            #str(side)
            if ( axes ) Axis(x, side=side)
 
            # This by() call coerces the piper object to a data.frame and then
            # the subset of the data.frame itself has to be coerced to a piper
            # object so that points() will call points.piper().
            if ( points ) {
              by(x, x@WaterType,
                 function(j, ...) {
                     points(x,j,...)
                 },
                 x, ... )
            }
            return(NULL)
          }
          )

# Coercion functions

as.data.frame.piper =
  function (x, ...)
  {
    as.data.frame.list( list( IDs=x@IDs, Ca=x@Ca, Mg=x@Mg, Cl=x@Cl,
                             SO4=x@SO4, group=x@group, WaterType=x@WaterType,
                             pt.pch=x@pt.pch, pt.col=x@pt.col ) )
  }


# TODO: work on this
setValidity("piper",
            function(object)
            {
              n <- sort(names(object))
              f <- c("Ca", "Cl", "Mg", "SO4")
              if ( identical( intersect(n,f) , f ) ) TRUE
              else paste("ERROR: missing items: ", setdiff(f,n), sep="")

              a <- length(object[[1]])
              nms <- names(object)
              x <- 0
              for (i in object){
                x <- x + 1
                nme <- nms[x]
              
                if (length(i) != object@Ca) {
                  stop("ERROR: lengths of items differ") }
                if (nme != "IDs" && nme != "WaterType") {
                #print(nme)
                  for (j in i) {
                    #print(paste("j=",j))
                    if ( ! is.numeric(j) || j < 1)
                      stop("ERROR: there is a non-numeric or negative number")
                  }
                }
              }

              over100 <- object$row.name[ ((object$Ca + object$Mg) > 100 |
                                           (object$Cl + object$SO4) > 100) ]
              if (length(over100) != 0) {
                print(paste("ERROR: row.name == ",over100,
                       ", Either cations or anions add up to more than 100%.",
                       sep=""))
                return(invisible())      # represses echo
              }
            }
            )

piper <- function(d, ...){
    new("piper", d, call=sys.call(), ...)
}

setMethod("show","piper",
          function(object){
              cat("Piper object: ")
              print(object@call)
              str(object)
          })
