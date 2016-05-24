# Author: Oscar Perpinan Lamigueiro oscar.perpinan@gmail.com
# Date :  October 2012
# Version 0.10
# Licence GPL v3

lineLabel <- function(line, label,
                      spar=.6, position = c('above', 'below'),
                      textloc = 'constantSlope',
                      col = add.text$col,
                      alpha = add.text$alpha,
                      cex = add.text$cex,
                      lineheight = add.text$lineheight,
                      font = add.text$font,
                      fontfamily = add.text$fontfamily,
                      fontface = add.text$fontface,
                      lty = add.line$lty,
                      lwd = add.line$lwd, 
                      col.line = add.line$col,
                      identifier = 'lineLabel',
                     ...){
  ## line is a list of Line objects
  ## label and line must be of equal length
  stopifnot(is.numeric(textloc) | is.character(textloc))
  stopifnot(length(line) == length(label))
  if (is.numeric(textloc)) stopifnot(length(line) == length(textloc))

  add.text <- trellis.par.get("add.text")
  add.line <- trellis.par.get("add.line")

  position <- match.arg(position)

  for (i in seq_along(line)) {
    subline <- line[[i]]
    lbl <- label[i]
    
    coords <- xy.coords(coordinates(subline))

    llines(coords,
           col = col.line, lty=lty, lwd=lwd,
           identifier = paste(identifier, "line", ## add ID?
             sep = "."))

    if (lbl !='') { ## Only execute if label is not empty
      ## Where the label is to be placed?
      ## Method "mixed" from lattice::panel.levelplot
      coords <- with(coords, data.frame(x, y))
      slopes <- diff(coords$y) / diff(coords$x)
      slopes[is.infinite(slopes)] <- NA
 
      rx <- range(coords$x)
      ry <- range(coords$y)

      stopifnot(is.character(textloc) | is.numeric(textloc))
      
      if (is.character(textloc)){
        textloc <- match.arg(textloc, c('minSlope', 'maxDepth', 'constantSlope'))

        idxLoc <- switch(textloc,
                          minSlope = which.min(abs(slopes)),
                          maxDepth = {
                            depth <- pmin(pmin(coords$x - rx[1], rx[2] - coords$x) / diff(rx), 
                                          pmin(coords$y - ry[1], ry[2] - coords$y) / diff(ry))
                            min(which.max(depth), length(slopes))
                            },
                          constantSlope = {
                            slopeSign <- sign(slopes)
                            zeros <- which(slopes==0)
                            for (i in zeros) {
                              i1 <- max(1, i-1)
                              i2 <- min(length(slopes), i + 1)
                              slopeSign[i] <- sign(sum(slopes[c(i1, i2)], na.rm=1))
                              }
                            rleSlopes <- rle(slopeSign)
                            idxCS <- max(which.max(rleSlopes$lengths) - 1, 1)
                            cumLengths <- cumsum(rleSlopes$lengths)
                            idxLim <- c(cumLengths[idxCS] +1, cumLengths[idxCS + 1])
                            idxLim[1]
                            ## idx1 <- cumsum(rleSlopes$lengths)[idxCS] + 1
                            ## idx2 <- idx1 + rleSlopes$lengths[idxCS + 1] - 1
                            ## which(coords[,1] == min(coords[c(idx1, idx2), 1]))
                            }
                          )
        } else {idxLoc <- textloc[i]}
      
      ## Label: number of characters, width and height
      gp <- gpar(col = col,
                 alpha = alpha,
                 cex = cex,
                 lineheight = lineheight,
                 fontfamily = fontfamily,
                 fontface = fontface)
      chars <- strsplit(lbl, '')[[1]]
      nchars <- length(chars)


      chWidth <- sapply(chars, function(ch){
        chGrob <- textGrob(ch, gp=gp)
        gWidth <- grobWidth(chGrob)
        convertWidth(gWidth, 'native', valueOnly=TRUE)
      })

      labLength <- sum(chWidth)

      chHeight <- sapply(chars, function(ch){
        chGrob <- textGrob(ch, gp=gp)
        gHeight <- grobHeight(chGrob)
        convertHeight(gHeight, 'native', valueOnly=TRUE)
      })
      labHeight <- max(chHeight)

      ##------------------------------------------------------------------##
      ## Starting at idxLoc, extract the segment of the line to be
      ## followed. It must be long enough for the text
      ##------------------------------------------------------------------##

      lengths <- LineLength(subline, sum=FALSE, longlat=FALSE)
      nParts <- length(lengths)
      lengthsLoc <- cumsum(lengths[idxLoc:nParts])

      if (labLength > lengthsLoc[nParts - idxLoc +1]) { ##The segment is not long enough for the label
        idxLoc <- 1
        segCoords <- apply(coords, 2, mean)
        slopes <- mean(slopes, na.rm=TRUE)
      } else {
        idxLength <- min(which(lengthsLoc >= labLength))
        segCoords <- coords[idxLoc:(idxLoc+idxLength),]
        segment <- Line(segCoords)
        lengthSeg <- LineLength(segment)
      }

      nX <- nrow(segCoords)
      if (!is.null(nX) && nX >= 4) { ## smooth.spline needs at least four distinct x values

        ##------------------------------------------------------------------##
        ## Resample the segment to obtain a collection of short segments.
        ##------------------------------------------------------------------##

        n <- floor(lengthSeg/min(chWidth)) * 1e3

        segSplineFun <- smooth.spline(segCoords, spar = spar) ## change spar to adjust the smooth level
        rxSeg <- range(segCoords$x)
        xsp <- seq(rxSeg[1], rxSeg[2], length=n) ## new x values to resample the segment
        segSpline <- data.frame(predict(segSplineFun, xsp))

        ## Each letter will be placed above a midpoint of this resampled segment
        midPoints <- 1/2 * (segSpline[1:(n-1), ] + segSpline[2:n, ])
        slopes <- diff(segSpline$y) / diff(segSpline$x)

        ##------------------------------------------------------------------##
        ## Parallel translation of the resampled segment. The distant is fixed
        ## by the height of the label.
        ##------------------------------------------------------------------##

        ## Distance between segments
        D <- labHeight
        ## Unitary perpendicular vector
        perp <- data.frame(x=-diff(segSpline$y), y=diff(segSpline$x))
        len <- apply(perp, 1, FUN=function(v)sqrt(sum(v^2)))
        perp <- perp/len

        ## Points above (+) or below (-) line

        if (position =='below') D <- abs(D) else D <- -abs(D)
        segOffset <- sweep(midPoints, 2, D * colMeans(perp))
        ##segOffset <- midPoints +  D*perp

        segSplineLine <- Line(segOffset) ##Line(segSpline)
        lengthSegSpline <- LineLength(segSplineLine, sum=FALSE)

        ## sp.lines(segSplineLine, col='blue')

        ##------------------------------------------------------------------##
        ## A character is placed above a collection of segments whose total
        ## length is higher than the character width. IDs is a vector of
        ## indices to subset the parallel segment.
        ##------------------------------------------------------------------##

        id <- IDs <- 1
        lengths <- lengthSegSpline
        nSegs <- length(lengths)
        for (i in 2:nchars){
          D <- cumsum(lengths)
          L <- 1/2*(chWidth[i-1] + chWidth[i])
          id <- min(which(D>=L))
          lengths <- lengths[(id+1):nSegs]
          IDs <- c(IDs, id)
        }
        IDs <- cumsum(IDs)

        pts <- SpatialPoints(segOffset)[IDs,] 
        coordLabs <- coordinates(pts)
        slopeLabs <- 180/pi * atan(slopes[IDs])

        tg <- textGrob(chars, coordLabs[,1], coordLabs[,2],
                       default.units='native',
                       rot=slopeLabs,
                       gp=gp)
        grid.draw(tg)

      } else { ## without 4 distinct x values
        coordLabs <- coords[idxLoc,]
        slopeLabs <- 180/pi * atan(slopes[idxLoc])
        just <- if (position=='below') 'top' else 'bottom'

        tg <- textGrob(lbl, coordLabs[,1], coordLabs[,2],
                       just=just,
                       default.units='native',
                       rot=slopeLabs,
                       gp=gp)
        grid.draw(tg)



      }
      ## lapply(1:nchars, FUN=function(i){
      ##   pg <- circleGrob(coordLabs[i,1], coordLabs[i,2], r=chWidth[i]/2,
      ##                    default.units='native', gp=gpar(col='red'))
      ##   grid.draw(pg)
      ##  })
    }}
}

setGeneric('sp.lineLabel', function(object, labels, byid=TRUE,...){standardGeneric('sp.lineLabel')})

setMethod('sp.lineLabel',
          signature=(object='Lines'),
          definition=function(object, labels, byid=TRUE, ...){
            if (missing(labels)) labels=object@ID
            if (isTRUE(byid)){
              lengths <- sapply(object@Lines, LineLength)
              idxMax <- which.max(lengths)
              labs <- rep('', length(object@Lines))
              labs[idxMax] <- labels
            } else {
              labs <- rep(labels, length(object@Lines))
            }
            lineLabel(object@Lines, label=labs, ...)
          }
          )

setMethod('sp.lineLabel',
          signature=(object='SpatialLines'),
          definition=function(object, labels, byid=TRUE,...){
            if (missing(labels)){
              ids <- sapply(object@lines, function(l)l@ID)
              labels <- ids
              names(labels) <- ids
            }
            for (i in seq_along(object)){
              line <- object@lines[[i]]
              idLabel <- match(line@ID, names(labels))
              lbl <- if (is.na(idLabel)) '' else as.graphicsAnnot(labels[idLabel])
              sp.lineLabel(line, lbl, byid=byid,...)
            }
          }
          )


setGeneric('label', function(object, text,...){standardGeneric('label')})

setMethod('label', signature=(object='SpatialLines'),
          definition=function(object, text,...){
            text <- as.character(text)
            nfeatures <- length(object)
            if (!identical(nfeatures, length(text))){
              stop('Number of labels must equal number of features.')
            }
            ids <- sapply(object@lines, function(l)l@ID)
            names(text) <- ids
            text
          })
