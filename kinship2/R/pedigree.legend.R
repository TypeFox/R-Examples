# Automatically generated from all.nw using noweb

pedigree.legend <- function (ped, labels = dimnames(ped$affected)[[2]],
    edges = 200, radius=NULL, location="bottomright", new=TRUE,
    density = c(-1, 35, 55, 25),  angle = c(90, 65, 40, 0), ...) 
{
   
    naff <- ncol(ped$affected)
    x <- rep(1,naff)
    
    # Defaults for plotting on separate page:
    ## start at the top, always counter-clockwise, black/white
    init.angle <- 90
    twopi <- 2 * pi
    col <- 1

    default.labels <- paste("affected-", 1:naff, sep='')
    if (is.null(labels)) labels <- default.labels
    
    ## assign labels to those w/ zero-length label
    whichNoLab <- which(nchar(labels) < 1)
    if(length(whichNoLab))
      labels[whichNoLab] <- paste("affected-", whichNoLab, sep='')

    
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    ## settings for plotting on a new page
    if(!new) {
      plot.new()
      
      pin <- par("pin")
      # radius, xylim, center, line-lengths set to defaults of pie()
      radius <- 1
      xlim <- ylim <- c(-1, 1)
      center <- c(0,0)
      llen <- 0.05
      
      if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L]/pin[2L]) * xlim
      else ylim <- (pin[2L]/pin[1L]) * ylim
      
      plot.window(xlim, ylim, "", asp = 1)
      
    } else {
      ## Settings to add to pedigree plot
      ## y-axis is flipped, so adjust angle and rotation
      init.angle <- -1*init.angle
      twopi <- -1*twopi

      ## track usr xy limits. With asp=1, it re-scales to have aspect ratio 1:1
      usr.orig <- par("usr")
      plot.window(xlim=usr.orig[1:2], ylim=usr.orig[3:4], "", asp=1)
      usr.asp1 <- par("usr")

      ## set line lengths
      llen <- radius*.15
     
      ## also decide on good center/radius if not given
      if(is.null(radius))
        radius <- .5
     
      ## get center of pie chart for coded
      pctusr <- .10*abs(diff(usr.asp1[3:4]))
      center = switch(location,
        "bottomright" = c(max(usr.asp1[1:2])-pctusr,max(usr.asp1[3:4])-pctusr),
        "topright" = c(max(usr.asp1[1:2])-pctusr, min(usr.asp1[3:4]) + pctusr),
        "bottomleft" =c(min(usr.asp1[1:2]) + pctusr, max(usr.asp1[3:4])-pctusr),
        "topleft" = c(min(usr.asp1[1:2]) + pctusr, min(usr.asp1[3:4]) + pctusr))
     
    }
    
    col <- rep(col, length.out = nx)
    border <- rep(1, length.out = nx)
    lty <- rep(1, length.out = nx)
    angle <- rep(angle, length.out = nx)
    density <- rep(density, length.out = nx)
  
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p), y = radius * sin(t2p))
    }
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        P$x <- P$x + center[1]
        P$y <- P$y + center[2]
        
        polygon(c(P$x, center[1]), c(P$y, center[2]), density = density[i],
                angle = angle[i], border = border[i], col = col[i],
                lty = lty[i])

        P <- t2xy(mean(x[i + 0:1]))
        if(new) {
          ## not centered at 0,0, so added center to x,y
          P$x <- P$x + center[1]
          P$y <- center[2] + ifelse(new, P$y, -1*P$y)
        }
        
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
          ## put lines
          lines(x=c(P$x, P$x + ifelse(P$x<center[1], -1*llen, llen)),
                y=c(P$y, P$y + ifelse(P$y<center[2], -1*llen, llen)))

          ##  put text just beyond line-length away from pie
          text(x=P$x + ifelse(P$x < center[1], -1.2*llen, 1.2*llen),
               y=P$y + ifelse(P$y < center[2], -1.2*llen, 1.2*llen),
               labels[i], xpd = TRUE, 
               adj = ifelse(P$x < center[1], 1, 0), ...)
        }
    }
    
    invisible(NULL)
}
