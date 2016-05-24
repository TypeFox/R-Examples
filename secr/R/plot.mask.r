#############################################################################
## package 'secr'
## plot.mask.R
## 2014-09-14 strip.legend() added
## 2014-09-20 strip.legend allows [,)
## 2014-09-30 covrange allws NA
## 2015-01-14 tweak to format legend text (breaks case)

## suggest: optionally suppress legend text
## suggest: optionally use axis()
## Help page for strip.legend
## add 'legend' argument to doc for plot.mask
#############################################################################

strip.legend <- function (xy, legend, col,
                         legendtype = c('breaks', 'intervals', 'other'),
                         tileborder = NA,
                         height = 0.5,
                         width = 0.06,
                         inset = 0.06,
                         text.offset = 0.02,
                         text.cex = 0.9,
                         xpd = TRUE,
                         scale = 1,
                         title = '',
                         box = NA,
                         box.col = par()$bg) {

  ## assumes equal-interval bars

  plotpixel <- function (yx, col) {
    ## plot a rectangle of color col centred at coordinates yx with half-sides dx,dy
    rect (yx[2]-dx, yx[1]-dy, yx[2]+dx, yx[1]+dy, col = col, density=-1, border = tileborder)
  }
  legendtype <- match.arg(legendtype)
  oldxpd <- par()$xpd
  on.exit(par(xpd = oldxpd))
  par(xpd = xpd)
  usr <- par()$usr
  mar <- par()$mar

  ##-------------------------------------
  ## tidy legend text
  legend <- gsub( '\\[', '', legend)
  legend <- gsub( '\\]', '', legend)
  legend <- gsub( '\\(', '', legend)
  legend <- gsub( '\\)', '', legend)
  breaks <- sapply(legend, strsplit, ',')
  breaks <- unique(unlist(breaks))
  if (all(!is.na(suppressWarnings(as.numeric(breaks)))))
      breaks <- sort(as.numeric(breaks)) * scale

  if (legendtype == 'breaks') {
    ## 2015-01-14
    ## legend <- as.character(breaks)
    legend <- format(breaks, trim = TRUE)
    ## legend <- sprintf(paste("%8.", dec, "f", sep=""), breaks)
    ## legend <- formatC( breaks, digits = 2, format = "g")
  }
  else if (legendtype == 'intervals') {
    legend <- paste(breaks[-length(breaks)], breaks[-1], sep = ' - ')
  }
  legwidth <- max(strwidth(legend))   ## required width in user units
  textoffset <- text.offset * (usr[2]-usr[1])
  ##-------------------------------------
  ## determine dimensions of each rectangle
  ncol <- length(col)
  wy <- (usr[4]-usr[3]) / ncol * height
  wx <- (usr[2]-usr[1]) * width
  dx <- wx / 2
  dy <- wy / 2
  ## --------------------------------------
  if (nchar(title) == 0)
    titleheight <- -wy/2
  else
    titleheight <- strheight(title, cex = text.cex)

  if (is.character(xy)) {
    ## assume "topright" placement
    if (!(xy %in% c("topleft", "topright", "bottomleft", "bottomright",
                    "left", "right")))
      stop ("only 'topleft', 'topright', 'bottomleft', 'bottomright',
             'left', and 'right' implemented")
    if (xpd & (par()$pin[1]>0) & ((mar[2] + mar[4])>0)){
      ## left and right margin width in user units
      userin <- (usr[2]-usr[1]) / par()$pin[1]  ## user units per inch
      xmarginwidth <- (par()$fin - par()$pin)[1]
      mx <- xmarginwidth * mar[c(2,4)]/(mar[2] + mar[4]) * userin
    }
    else {
      mx <- c(0,0)
    }
    ## programmed locations
    insetx <- inset * (usr[2]-usr[1])
    insety <- inset * (usr[4]-usr[3])
    if (xy == 'topright')
      xy <- c(usr[2] - legwidth - wx - textoffset + mx[2] - insetx, usr[4]-insety)
    else if (xy == 'topleft')
      xy <- c(usr[1] - mx[1] + insetx, usr[4]-insety)
    else if (xy == 'bottomright')
      xy <- c(usr[2] - legwidth - wx - textoffset + mx[2] - insetx, usr[3]+insety+wy*ncol)
    else if (xy == 'bottomleft')
      xy <- c(usr[1] - mx[1] + insetx, usr[3]+insety+wy*ncol)
    else if (xy == 'left')
      xy <- c(usr[1] - mx[1] + insetx,
              (usr[3]+usr[4])/2 + wy*ncol/2)
    else if (xy == 'right')
      xy <- c(usr[2] - legwidth - wx - textoffset + mx[2] - insetx,
              (usr[3]+usr[4])/2 + wy*ncol/2)

  }
  else {
    xy <- unlist(xy)
    xy <- xy + c(textoffset, - titleheight - textoffset)
  }

  ## locate vertical centres and bounds

  centres <- seq(xy[2]-wy/2, by = -wy, length.out = ncol)
  boundsy <- seq(xy[2], by = -wy, length.out = ncol+1)
  centres <- rev(centres)   ## low to high
  boundsy <- rev(boundsy)   ## low to high

  ## optional frame
  if (!is.na(box))
    if (is.logical(box)) box <- if(box) par()$fg else NA
  if (!is.na(box)) {
    rect(xy[1]-textoffset,
         centres[1] - wy/2 - textoffset ,
         xy[1]+ wx + 2*textoffset + legwidth,
         centres[ncol] + titleheight + textoffset + wy,
         border = box,
         col = box.col
    )
  }

  ## strip
  yx <- lapply(centres, c, xy[1]+wx/2)
  tmp <- mapply(plotpixel, yx, col)

  ## add legend header
  text(xy[1] + (wx + textoffset + legwidth)/2,
       centres[ncol] + strheight(title, cex = text.cex) + wy,
       title, adj = 0.5, cex = text.cex )

  ## add legend text to right of this rectangle
  alt <- rep(c(TRUE,FALSE), length.out = length(legend))
  if (legendtype == 'other')
    text (rep(xy[1] + wx + textoffset, ncol), centres, legend,
        adj = 0, cex = text.cex)
  else if (sum(strheight(legend, cex = text.cex)) < diff(range(centres))) {
      if (legendtype == 'breaks')
          text (rep(xy[1] + wx + textoffset, ncol+1), boundsy, legend,
                adj = 0, cex = text.cex)
      else if (legendtype == 'intervals')
          text (rep(xy[1] + wx + textoffset, ncol), centres, legend,
                adj = 0, cex = text.cex)
  }
  else if (sum(strheight(legend[alt], cex = text.cex)) < diff(range(centres))) {
      if (legendtype == 'breaks')
          text (rep(xy[1] + wx + textoffset, ncol+1)[alt], boundsy[alt], legend[alt],
                adj = 0, cex = text.cex)
      else if (legendtype %in% c('intervals', 'other'))
          text (rep(xy[1] + wx + textoffset, ncol)[alt], centres[alt], legend[alt],
                adj = 0, cex = text.cex)
  }
  else {
      if (legendtype == 'breaks')
          text (rep(xy[1] + wx + textoffset, 2), boundsy[c(1,ncol+1)],
                legend[c(1,ncol+1)], adj = 0, cex = text.cex)
      else if (legendtype %in% c('intervals', 'other'))
          text (rep(xy[1] + wx + textoffset, 2), centres[c(1,ncol)],
                legend[c(1,ncol)], adj = 0, cex = text.cex)
  }

  invisible(c(xy[1], xy[1]+wx, xy[2], xy[2] - ncol * wy))

}
#############################################################################

plot.mask <- function(x, border = 20, add = FALSE, covariate = NULL,
      axes = FALSE, dots = TRUE, col = 'grey', breaks = 10, meshcol = NA,
      ppoly = TRUE, polycol = 'red', legend = TRUE, ...)
{
    if (ms(x)) {
        ## 2013-02-12 pass all arguments
        lapply (x, plot.mask, border = border, add = add, covariate = covariate,
                axes = axes, dots = dots, col = col, breaks = breaks, meshcol =
                meshcol, ppoly = ppoly, polycol = polycol, legend = legend, ...)
    }
    else {
        buff <- c(-border,+border)
        if (!add) {
            eqscplot (x$x, x$y,
                      xlim = range(x$x) + buff,
                      ylim = range(x$y) + buff,
                      xlab = '', ylab = '',
                      axes = axes, type = 'n')
        }

        if (!is.null(attr(x,'polygon')) & ppoly) {
            poly <- attr(x,'polygon')
            if (class(poly) == "SpatialPolygonsDataFrame") {
# plot(poly, col = polycol, add = TRUE)
# poor control of colours
                plot(poly, add = TRUE)
            }
            else
                polygon (poly, col = polycol, density = 0)
        }

        if (is.null(covariate))
            covfactor <- factor(1)
        else {
            if (is.factor(covariates(x)[,covariate]))
                covfactor <- covariates(x)[,covariate]
            else {
              covvalue <- covariates(x)[,covariate]
              if (length(breaks) == 1) {
                ## covrange <- range(covvalue)
                covrange <- range(covvalue, na.rm = TRUE)   ## na.rm 2014-09-30
                rough <- seq(covrange[1], covrange[2], length.out = breaks+1)
                breaks <- pretty(rough, n = breaks)
              }
              ## include.lowest = TRUE added 2014-09-20
              covfactor <- cut ( covvalue, breaks = breaks, include.lowest = TRUE)
            }
        }
        ncolour <- length(levels(covfactor))
        if (length(col) < ncolour) {
            ## col <- heat.colors(ncolour)   # old default set
            if (length(col) > 1)
                warning ("too few colours; using terrain.colors(", ncolour, ")")
            col <- terrain.colors(ncolour)   # new default set 2.9.0
        }
        cols <- col[as.numeric(covfactor)]
        allargs <- list(...)
        if (dots) {
                args <- list(x= as.data.frame(x), pch = 16, cex = 0.8
                                , col = cols, type = 'p')
                dotsargs <- allargs[names(allargs) %in% c('pch','cex','type')]
                args <- replace(args, names(dotsargs), dotsargs)
                do.call(points, args)
        }
        else {
            pixelsize <- attr(x,'spacing')
            ## dx <- c(-0.5, -0.5, +0.5, +0.5) * pixelsize
            ## dy <- c(-0.5, +0.5, +0.5, -0.5) * pixelsize
            dx <- pixelsize / 2
            dy <- pixelsize / 2
            plotpixel <- function (xy) {
                rect (xy[1]-dx, xy[2]-dy, xy[1]+dx, xy[2]+dy, col = col[xy[3]],
                      density=-1, border = meshcol)
            }
            apply(cbind(x,as.numeric(covfactor)),1,plotpixel)
        }
        if (legend & !is.null(covariate)) {
            legendtext <- levels(covfactor)[1:ncolour]
            if (dots) {
                args <- formals(legend)
                newargs <- list(x = 'right', legend = rev(legendtext), pch = 16,
                       col = rev(col[1:ncolour]), title = covariate)
                args <- replace(args, names(newargs), newargs)
                args <- replace(args, names(dotsargs), dotsargs)           
                if ('xy' %in% names(allargs))
                    args$x <- allargs$xy
                do.call('legend', args)
            }
            else { 
                args <- formals(strip.legend)
                newargs <- list(xy = 'right', col = col[1:ncolour],
                                legend = legendtext, tileborder = meshcol,
                                title = covariate)
                if (is.factor(covariates(x)[,covariate])) {
                  newargs$legendtype <- 'other'
                  newargs$height <- min(1, length(legendtext) * 0.06)
                }
                args <- replace(args, names(newargs), newargs)
                args <- replace(args, names(allargs), allargs)
                do.call(strip.legend, args)
            }
        }

        if (!is.null(covariate))
            invisible(levels(covfactor)[1:ncolour])
    }
}
###############################################################################
