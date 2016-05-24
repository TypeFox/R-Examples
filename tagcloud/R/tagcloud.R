tagcloud.debug <- FALSE
#tagcloud.debug <- TRUE
debugprf  <- function( ... ) if( tagcloud.debug ) print( sprintf( ... ))
debugpr   <- function( ... ) if( tagcloud.debug ) print( ... )
debugcatf <- function( ... ) if( tagcloud.debug ) cat( sprintf( ... ))

# calculate cex from floor, ceiling and value vector
.calc.cex <- function( x, floor, ceiling, wmin= NULL, wmax= NULL ) {

  if( is.null( wmin ) ) wmin <- min( x )
  if( is.null( wmax ) ) wmax <- max( x )

  ret <- (floor + (ceiling - floor) * (x - wmin) / wmax) 
  ret
}

.calc.sizes <- function( tags, boxes, family ) {

  if( length( family ) < 2 ) family <- rep( family, length( tags ) )

  for( i in 1:length( tags ) ) {
    boxes[ i, "w" ] <- strwidth( paste0( "X", tags[i] ), cex= boxes[i, "cex"], family= family[ i ] )
    boxes[ i, "h" ] <- 1.35 * strheight( tags[i], cex= boxes[i, "cex"], family= family[ i ] )
  }
  boxes[ , "s" ] <- boxes[ , "w" ] * boxes[ , "h" ]

  boxes
}


.auto.scale <- function( boxes ) {

  usr <- par( "usr" )
  plot.surface <- ( usr[2] - usr[1] ) * ( usr[4] - usr[3] )
  box.surface  <- sum( boxes[,"s"] )
  scale <- sqrt( 0.4 * plot.surface / box.surface )
  if( nrow( boxes ) < 20 ) scale <- scale * 0.8
  debugprf( "box=%.2f plot=%.2f ratio=%.2f scale=%.2f",
    box.surface, plot.surface, plot.surface / box.surface, scale )
  return( scale )
}

# reduce space between non-overlapping boxes
.squeeze.boxes <- function( boxes ) {

  c <- .box.center( boxes )
  n <- nrow( boxes )

  if( tagcloud.debug ) {
    debugprf( "before: bb: %.2f, %.2f -> %.2f, %.2f", 
      min( boxes[,"x"] ), min( boxes[,"y"] ), max( boxes[,"x"] + boxes[,"w"] ), max( boxes[,"y"] + boxes[,"h"] ) )
    debugprf( "surface index before squeezing: %.3f", boxes.ratio( boxes ) )
  }

  for( nn in 1:25 ) {
    ydist <- boxes[,"y"] + boxes[,"h"] / 2 - c[2]
    ord.y <- order( abs( ydist ) )

    for( o in 2:length( ord.y ) ) {
      b <- ord.y[o]
      y.step <- ydist[b] / 15 
      
      for( i in 1:20 ) {
        if( is.overlap( c( boxes[b,"x"], boxes[b,"y"] - y.step, boxes[b,"w"], boxes[b,"h"] ), boxes[-b,,drop=F] ) ) break ;
        boxes[b,"y"] <- boxes[b,"y"] - y.step
      }
    }
  }

  for( nn in 1:25 ) {
    xdist <- boxes[,"x"] + boxes[,"w"] / 2 - c[1]
    ord.x <- order( abs( xdist ) )
    for( o in 2:length( ord.x ) ) {
      b <- ord.x[o]
      x.step <- xdist[b] / 150
      for( i in 1:20 ) {
        if( is.overlap( c( boxes[b,"x"] - x.step, boxes[b,"y"], boxes[b,"w"], boxes[b,"h"] ), boxes[-b,,drop=F] ) ) break ;
        boxes[b,"x"] <- boxes[b,"x"] - x.step
      }
    }
  }


  if( tagcloud.debug ) {
    debugprf( "surface index before squeezing: %.3f", boxes.ratio( boxes ) )
    debugprf( "after: bb: %.2f, %.2f -> %.2f, %.2f", 
      min( boxes[,"x"] ), min( boxes[,"y"] ), max( boxes[,"x"] + boxes[,"w"] ), max( boxes[,"y"] + boxes[,"h"] ) )
  }

  return( boxes )
}

.box.center <- function( boxes ) {

  x.min <- min( boxes[,"x" ] )
  y.min <- min( boxes[,"y" ] )
  x.max <- max( boxes[,"x" ] + boxes[,"w"])
  y.max <- max( boxes[,"y" ] + boxes[,"h"])
  
  xc <- x.min + ( x.max - x.min ) / 2 
  yc <- y.min + ( y.max - y.min ) / 2 

  return( c( xc, yc ) )

}

.center.boxes <- function( boxes ) {
  usr <- par( "usr" )

  x0 <- usr[1] + ( usr[2] - usr[1] ) / 2
  y0 <- usr[3] + ( usr[4] - usr[3] ) / 2

  c <- .box.center( boxes )

  boxes[,"x"] <- boxes[,"x"] - c[1] + x0
  boxes[,"y"] <- boxes[,"y"] - c[2] + y0
  return( boxes )
}

.get.min.scale <- function( boxes ) {
  
  tot.w <- max( boxes[,"x"] + boxes[,"w"] ) - min( boxes[,"x"] ) 
  tot.h <- max( boxes[,"y"] + boxes[,"h"] ) - min( boxes[,"y"] ) 

  usr <- par( "usr" )

  scale.w <- ( usr[2] - usr[1] ) / tot.w
  scale.h <- ( usr[4] - usr[3] ) / tot.h

  debugprf( "tot.w= %.2f, usr.w= %.2f, scale.w=%.2f", tot.w, usr[2] - usr[1], scale.w )
  debugprf( "tot.h= %.2f, usr.h= %.2f, scale.h=%.2f", tot.h, usr[4] - usr[3], scale.h )

  return( min( scale.w, scale.h ) )

}

.fit.boxes <- function( tags, boxes, vertical, family ) {
  stop()

  debugprf( "*********************************" )
  debugprf( "initial surface index: %.3f", boxes.ratio( boxes ) )
  c <- .box.center( boxes )
  scale <- .get.min.scale( boxes )
  debugprf( "scale=%.2f", scale )
  if( any.overlap( boxes ) ) stop( "boxes overlap" )

  max.iter <- 20
  f1 <- 1
  f2 <- 0
  f  <- 1
  scale.new <- -99
  if( scale < 1 ) stop( "not yet" )
  while( max.iter > 0 ) {

    debugprf( "iter=%d f1=%.2f f2=%.2f f=%.2f scale=%.2f", max.iter, f1, f2, f, scale.new )

    if( f2 == 0 ) f <- f * 2
    else          f <- ( f2 + f1 ) / 2
    debugprf( "        f1=%.2f f2=%.2f f=%.2f scale=%.2f", f1, f2, f, scale.new )

    boxes.tmp <- boxes
    boxes.tmp[, "cex"] <- boxes.tmp[,"cex"] * f
    boxes.tmp <- .calc.sizes( tags, boxes.tmp, family )
    boxes.tmp[ vertical, c( "w", "h" ) ] <- boxes.tmp[ vertical, c( "h", "w" ) ]
    scale.new <- max( c( boxes.tmp[,"w"] / boxes[,"w"], boxes.tmp[,"h"] / boxes[,"h"] ) )
    boxes.tmp[,"x"]   <- c[1] + ( boxes.tmp[,"x"] - c[1] ) * scale.new
    boxes.tmp[,"y"]   <- c[2] + ( boxes.tmp[,"y"] - c[2] ) * scale.new

    if( scale.new > scale ) {
      debugprf( "** over scale" )
      f2 <- f
    } else {
      debugprf( "** under scale" )
      f1 <- f
    }


    max.iter <- max.iter - 1 
  }

  boxes.orig <- boxes
  boxes[, "cex"] <- boxes[,"cex"] * f1
  boxes <- .calc.sizes( tags, boxes, family )
  boxes[ vertical, c( "w", "h" ) ] <- boxes[ vertical, c( "h", "w" ) ]
  scale.new <- max( c( boxes[,"w"] / boxes.orig[,"w"], boxes[,"h"] / boxes.orig[,"h"] ) )
  boxes[,"x"]   <- c[1] + ( boxes[,"x"] - c[1] ) * scale.new
  boxes[,"y"]   <- c[2] + ( boxes[,"y"] - c[2] ) * scale.new

  return( boxes )
}


boxes.width  <- function( boxes ) max( boxes[,"x"] + boxes[,"w"] ) - min( boxes[,"x"] ) 
boxes.height <- function( boxes ) max( boxes[,"y"] + boxes[,"h"] ) - min( boxes[,"y"] ) 

boxes.ratio <- function( boxes ) {

  tot.surface <- boxes.width( boxes ) * boxes.height( boxes )
  tag.surface <- sum( boxes[,"w"] * boxes[,"h"] )

  tag.surface / tot.surface
}

auto.wh <- function( boxes, mult= 2 ) {
  #ds <- dev.size()
  ds <- par( "pin" )
  aspect <- ds[2] / ds[1]
  w <- sum( boxes[, "w" ] ) / mult
  h <- sum( boxes[, "h" ] ) / mult
  if( h > w ) w <- h / aspect 
  else        h <- w * aspect 
  ret <- c( w, h, aspect, ds )
  names( ret ) <- c( "w", "h", "aspect", "pinw", "pinh" )
  return( ret )
}


# fit using normal distribution
algorithm.normal <- function( boxes ) {

  i <- 20 
  while( i > 0 ) {
    boxes[,"x"] <- rnorm( nrow( boxes ), sd= 0.1 ) - boxes[,"w"]/2 
    boxes[,"y"] <- rnorm( nrow( boxes ), sd= 150 )
    i <- i - 1
    if( ! any.overlap( boxes ) ) {
      debugprf( "not bad!" )
      break ;
    }
  }

  boxes <- .squeeze.boxes( boxes )
  boxes <- .center.boxes( boxes )

  return( boxes )
}


# attempt a uniform distribution 
algorithm.random <- function( boxes ) {

  wh <- auto.wh( boxes )
  debugpr( wh )
  w <- wh[ "w" ]
  h <- wh[ "h" ]
  # this is our target aspect
  aspect <- wh[ "aspect" ]

  debugprf( "w=%.2f, h=%.2f, aspect= %.2f", w, h, wh[ "aspect" ] )

  boxes.new <- NULL

  for( i in 1:nrow( boxes ) ) {
    # how many tries
    n <- 100

    overlap <- TRUE
    while( overlap & n > 0 ) {
      x <- runif( 1, min= 0, max= w - boxes[i,"w"] )
      y <- runif( 1, min= 0, max= h - boxes[i,"h"] )
      if( i == 1 ) {
        overlap <- FALSE
        break 
      }

      if( ! any.overlap( rbind( boxes.new[,c( "x", "y", "w", "h")], c( x, y, boxes[i, "w"], boxes[i, "h" ] ) ) ) ) 
        overlap <- FALSE

      n <- n - 1
    }

    if( ! overlap ) {
      boxes[i,"x"] <- x
      boxes[i,"y"] <- y
      boxes.new <- rbind( boxes.new, boxes[i,c("x","y","w","h")] )
    } else {
      debugprf( "problems with word %d", i )
      stop( "Failed attempt! Try again" )
    }

  }

  boxes <- .squeeze.boxes( boxes )

  debugprf( "surface index before optimization: %.3f", boxes.ratio( boxes ) )
  debugprf( "total width before optimization: %.2f", boxes.width( boxes ) )
  debugprf( "total height before optimization: %.2f", boxes.height( boxes ) )
  debugprf( "current aspect: %.2f", boxes.height( boxes ) / boxes.width( boxes ) )

  for( loop in 1:20 ) {

    a <- boxes.height( boxes ) / boxes.width( boxes )
    debugprf( "loop= %d aspect= %.2f", loop, a )
    c <- .box.center( boxes )

    if( a < aspect ) {
      i <- which.max( sapply( 1:nrow( boxes ), function( x ) max( c[1] - boxes[x,"x"], c[1] - ( boxes[x,"x"] + boxes[x,"w"] ) ) ) )
    } else {
      i <- which.max( sapply( 1:nrow( boxes ), function( x ) max( c[2] - boxes[x,"y"], c[2] - ( boxes[x,"y"] + boxes[x,"h"] ) ) ) )
    }

    debugprf( "moving %d", i )
    c.2 <- .box.center( boxes[-i,] )
    w.2 <- boxes.width( boxes[-i,] )
    h.2 <- boxes.height( boxes[-i,] )
    a.2 <- h.2 / w.2
    debugprf( "remaining aspect: %.2f", a.2 )

    # make sure width is at least width of the i box, plus something
    if( w.2 < 1.1 * boxes[i, "w"] ) w.2 <- 1.1 * boxes[i, "w" ]
    if( h.2 < 1.1 * boxes[i, "h"] ) h.2 <- 1.1 * boxes[i, "h" ]

    # expand the area to get the correct aspect
    if( a.2 < aspect ) h.2 <- w.2 * aspect
    else               w.2 <- h.2 / aspect

    for( attempt in 1:5 ) {
      x <- runif( 1, min= c.2[1] - w.2 / 2, max= c.2[1] + w.2 / 2 )
      y <- runif( 1, min= c.2[2] - h.2 / 2, max= c.2[2] + h.2 / 2 )
      box.tmp <- rbind( boxes[-i, c( "x", "y", "w", "h" )], c( x, y, boxes[i, "w"], boxes[i, "h" ] ) )
      if( ! any.overlap( box.tmp ) ) {
        debugprf( "hah! (%d)", attempt )
        boxes[i,"x"] <- x
        boxes[i,"y"] <- y
        break ;
      }
    }
    boxes <- .squeeze.boxes( boxes )
  }

  debugprf( "surface index after optimization: %.3f", boxes.ratio( boxes ) )
  debugprf( "total width after optimization: %.2f", boxes.width( boxes ) )
  debugprf( "total height after optimization: %.2f", boxes.height( boxes ) )
  debugprf( "aspect after optimization: %.2f", boxes.height( boxes ) / boxes.width( boxes ) )

  #boxes <- .squeeze.boxes( boxes )
  boxes <- .center.boxes( boxes )
  return( boxes )
}


algorithm.list <- function( boxes, meta, centered= F ) {

  usr   <- par( "usr" )
  tot.h <- usr[4]-usr[3]
  tot.w <- usr[2]-usr[1]
  meta$vertical <- FALSE
  boxes <- .calc.sizes( meta$tags, boxes, meta$family )
  m <- 1.05

  while( sum( boxes[,"h"] * m ) < tot.h && max( boxes[,"w"] ) < tot.w ) {
    boxes[,"cex"] <- boxes[,"cex"] * 1.1
    boxes <- .calc.sizes( meta$tags, boxes, meta$family )
  }


  while( sum( boxes[,"h"] * m ) > tot.h || max( boxes[,"w"] ) > tot.w ) {
    boxes[,"cex"] <- boxes[,"cex"]*0.9
    boxes <- .calc.sizes( meta$tags, boxes, meta$family )
  }


  boxes[1,"y"] <- 0

  for( i in 2:nrow( boxes ) ) {
    boxes[i,"y"] <- boxes[i-1,"y"] - m * boxes[i,"h"]
  }

  if( centered ) {
    boxes[,"x"] <- 0 - boxes[,"w"] / 2 
  } else {
    boxes[,"x"] <- 0
  }


  boxes <- .center.boxes( boxes )
  return( boxes )
}


algorithm.ulam <- function( boxes ) {

  debugprf( "*** algorithm ulam" )

  wh <- auto.wh( boxes )
  w  <- wh[ "w" ]
  h  <- wh[ "h" ]
  aspect <- wh[ "aspect" ]

  max.r <- 10 * sqrt( w^2 + h^2 ) 
  debugprf( "w=%.2f, h=%.2f, aspect=%.2f, max.r=%.2f", w, h, aspect, max.r )

  r.step <- min( boxes[,"w"] ) / 15

  boxes[1,"x"] <- 0 - boxes[1,"w"]/2
  boxes[1,"y"] <- 0 - boxes[1,"h"]/2
  boxes.new <- boxes[1,,drop=F]

  for( i in 2:nrow( boxes ) ) {
    debugcatf( "\rCalculating, %d %% done", as.integer( 100 * i / nrow( boxes ) ) )
    r <- r.step
    x <- 0 - boxes[i,"w"]/2
    y <- 0 - boxes[i,"h"]/2

    # random initial direction
    d.r <- 0
    dir <- sample( 1:8, 1 )
    dir <- c( 0, 1, -1, 0, 0, -1, 1, 0,
              0,-1,  1, 0, 0,  1,-1, 0 )[(dir*2-1):(dir*2)]

    overlap <- TRUE

    params <- list( x=x, y=y, w=boxes[i,"w"], h=boxes[i,"h"],
                    dir1=dir[1], dir2=dir[2], rstep= r.step,
                    aspect= aspect, maxr= max.r, max.iter= 2000000 ) 

    test <- run.ulam( params, boxes[1:(i-1),c("x","y","w","h"),drop=F] ) 

    boxes[i,"x"] <- test[1]
    boxes[i,"y"] <- test[2]
    if( any.overlap( boxes[1:i,] ) ) stop() ;
  }

  debugcatf( "\n" )
  boxes <- .center.boxes( boxes )
  return( boxes )
}



# spiral algorithm, like in wordle
algorithm.spiral <- function( boxes ) {

  debugprf( "*** algorithm spiral" )
  wh <- auto.wh( boxes )
  w  <- wh[ "w" ]
  h  <- wh[ "h" ]
  aspect <- wh[ "aspect" ]

  max.r <- 10 * sqrt( w^2 + h^2 ) 
  debugprf( "w=%.2f, h=%.2f, aspect=%.2f, max.r=%.2f", w, h, aspect, max.r )

  boxes.new <- boxes[1,]

  r.step <- min( boxes[,"h"] ) / 15
  a.step <- pi / 90

  boxes[1,"x"] <- 0 - boxes[1,"w"]/2
  boxes[1,"y"] <- 0 - boxes[1,"h"]/2
  boxes.new <- boxes[1,,drop=F]

  dir <- 1

  for( i in 2:nrow( boxes ) ) {
    if( runif( 1 ) < 0.5 ) dir <- dir * -1 
    debugcatf( "\rCalculating, %d %% done", as.integer( 100 * i / nrow( boxes ) ) )
    angle  <- runif( 1, 0, 2 * pi )
    r <- mean( boxes[,"w"] / 3 )
    overlap <- TRUE
    x <- 0 - boxes[i,"w"]/2
    y <- 0 - boxes[i,"h"]/2

    params <- list( x=x, y=y, w=boxes[i,"w"], h=boxes[i,"h"], r= r,
                    dir=dir, rstep= r.step, angle= angle, astep= a.step,
                    aspect= aspect, maxr= max.r, max.iter= 2000000 ) 

    test <- run.spiral( params, boxes[1:(i-1),c("x","y","w","h"),drop=F] ) 

    boxes[i,"x"] <- test[1]
    boxes[i,"y"] <- test[2]
    if( any.overlap( boxes[1:i,] ) ) {
      printf( "overlap detected, i=%d", i )
      print( boxes[1:i,] ) 
      stop() 
    }
  }

  debugcatf( "\n" )
  boxes <- .center.boxes( boxes )
  return( boxes )
}

.boxesbox <- function( boxes ) {

  return( c( x0= min( boxes[,"x"] ),
             y0= min( boxes[,"y"] ),
             x1= max( boxes[,"x"] + boxes[,"w"] ),
             y1= max( boxes[,"y"] + boxes[,"h"] ) ) )

}

algorithm.snake <- function( boxes, meta ) {
  debugprf( "*** algorithm snake" )

  x0 <- 0
  y0 <- 0
  dirs <- c( "d", "l", "u", "r" )
  # srt  <- c(  d=270, l=180, u=90, r=0 )
  srt  <- c(  d=270, l=0, u=90, r=0 )

  meta$vertical <- FALSE
  boxes <- .calc.sizes( meta$tags, boxes, meta$family )
  
  boxes[1,"x"] <- boxes[1,"y"] <- 0
  n <- 1
  d <- 'd'
  px <- boxes[1,"w"]
  py <- boxes[1,"h"]

  prev.block <- c( x0=0, y0=0, x1=boxes[1,"w"], y2=boxes[1,"h"] )

  for( i in 2:nrow( boxes ) ) {

    if( d == 'd' ) {
      if( py - prev.block["y0"] < boxes[i,"w"] / 2 ) {
        d <- 'l'
        prev.block <- .boxesbox( boxes[1:(i-1),,drop=F] )
        px <- prev.block["x1"]
        py <- prev.block["y0"]
      }
    } else if( d == 'l' ) {
      if( px - prev.block["x0"] < boxes[i,"w"] / 2 ) {
        d <- 'u'
        prev.block <- .boxesbox( boxes[1:(i-1),,drop=F] )
        px <- prev.block["x0"]
        py <- prev.block["y0"]
      }
    } else  if( d == 'u' ) {
      if( prev.block["y1"] - py < boxes[i,"w"] / 2 ) {
        d <- 'r'
        prev.block <- .boxesbox( boxes[1:(i-1),,drop=F] )
        px <- prev.block["x0"]
        py <- prev.block["y1"]
      }
    } else if( d == 'r' ) {
      if( prev.block["x1"] - px < boxes[i,"w"] / 2 ) {
        d <- 'd'
        prev.block <- .boxesbox( boxes[1:(i-1),,drop=F] )
        px <- prev.block["x1"]
        py <- prev.block["y1"]
      }
    }

    if( d %in% c( 'd', 'u' ) ) boxes[i,c("h","w")] <- boxes[i,c("w","h")]
    boxes[i,"srt"] <- srt[d]

    if( d == 'd' ) {
      boxes[i,"x"] <- px
      boxes[i,"y"] <- py - boxes[i,"h"]

      px <- px
      py <- boxes[i,"y"]
    } else if( d == 'l' ) {
      boxes[i,"x"] <- px - boxes[i,"w"]
      boxes[i,"y"] <- py - boxes[i,"h"]

      px <- boxes[i,"x"]
      py <- py
    } else if( d == 'u' ) {
      boxes[i,"x"] <- px - boxes[i,"w"]
      boxes[i,"y"] <- py

      px <- px
      py <- py + boxes[i,"h"]
    } else if( d == 'r' ) {
      boxes[i,"x"] <- px
      boxes[i,"y"] <- py

      px <- px + boxes[i,"w"]
      py <- py
    }


    if( i == nrow( boxes ) ) break ;

    if( d == 'd' ) {
      if( py < prev.block["y0"] ) {
        py <- prev.block["y0"]
        prev.block <- .boxesbox( boxes[1:i,,drop=F] )
        d <- 'l'
      }
    } else if( d == 'l' ) {
      if( px < prev.block["x0"] ) {
        px <- prev.block["x0"]
        prev.block <- .boxesbox( boxes[1:i,,drop=F] )
        d <- 'u'
      }
    } else if( d == 'u' ) {
      if( py > prev.block["y1"] ) {
        py <- prev.block["y1"]
        prev.block <- .boxesbox( boxes[1:i,,drop=F] )
        d <- 'r'
      }
    } else if( d == 'r' ) {
      if( px > prev.block["x1"] ) {
        px <- prev.block["x1"]
        prev.block <- .boxesbox( boxes[1:i,,drop=F] )
        d <- 'd'
      }
    }
  }

  boxes <- .center.boxes( boxes )
  return( boxes )
}



# create a usable tagcloud object
.tagcloud.new.object <- function( boxes, meta, algorithm, scale ) {

  boxes <- as.data.frame( boxes )
  boxes <- cbind( meta, boxes )
  class( boxes ) <- c( "tagcloud", class( boxes ) )
  attr( boxes, "algorithm" ) <- algorithm
  attr( boxes, "scale" ) <- scale
  return( boxes )
}

# ----------------------------------------------------------------------
# main tagcloud functions
# ----------------------------------------------------------------------


#' Tag and Word Clouds
#' 
#' Functions to create and display plots called tag clouds, word clouds or
#' weighted lists, in which a usually large number of words is displayed in
#' size correlated with a numerical value (for example, frequency in a text or
#' enrichment of a GO term). This makes it easier to visualize the prominence
#' of certain tags or words. Also, it looks nice.
#' 
#' 
#' The package \code{tagcloud} creates and plots tag clouds (word clouds).  The
#' algorithms in the package have been designed specifically with long tags
#' (such as GO Term descriptions) in mind.
#' 
#' \subsection{Term ordering}{
#' The available arguments are as follows:
#' 
#' \itemize{
#'    \item size -- tags are ordered by size, that is, their effective
#'    width multiplied by their effective height. Default.  
#'    \item keep -- keep the order from the list of words provided 
#'    \item random -- randomize the tag list
#'    \item width -- order by effective screen width 
#'    \item height -- order by effective screen height 
#' }
#' 
#' By default, prior to plotting terms are ordered by size.
#' 
#' }
#' 
#' \subsection{Algorithms}{ 
#' There are four algorithms for placing tags on the
#' plot implemented in tagcloud.
#' 
#' \itemize{ 
#' \item oval -- creates an oval cloud.  
#' \item fill -- an attempt will
#' be made to fill the available space 
#' \item random -- randomly distribute tags
#' over the available space.  This algorithm is slow and not very effective
#' \item snake -- tags are placed clockwise around the first tag to plot
#' \item list -- create a list, one tag directly beneath another, justified left
#' \item clist -- create a list, one tag directly beneath another, centered 
#' }
#' 
#' Algorithms \code{oval}, \code{fill} and \code{random} attempt to fill the
#' available space by adjusting the scaling factor for the font sizes. }
#' 
#' \subsection{Calculation of tag sizes}{
#' Placing tags such that the empty space
#' between the tags is minimized poses a non-trivial problem, because the
#' effective bounding box of a displayed text is not linearly dependent on the
#' \code{cex} parameter.
#' 
#' In tagcloud, first a \code{cex} parameter is calculated for each tag
#' separately, based on the parameters \code{floor}, \code{ceiling} and the
#' vector of weights. Note that all weights smaller than \code{wmin} are
#' replaced by \code{wmin} and all weights larger than \code{wmax} are replaced
#' by \code{wmax}. Then, effective heights and widths of the tags to be
#' displayed are calculated using the \code{\link{strwidth}} and
#' \code{\link{strheight}} functions.
#' 
#' Unless the argument \code{scale} is different from "auto", a scaling
#' parameter for \code{cex} is automatically calculated based on the current
#' area of the tags and the available plotting area. This usually results in a
#' reasonable plot, but is neither guaranteed to occupy all of the available
#' space without margins, nor that no tag will cross the view port.
#' 
#' }
#' 
#' @aliases tagcloud tagcloud-package tagcloud-class plot.tagcloud
#' summary.tagcloud print.tagcloudsummary print.tagcloud
#' @param tags A character vector containing words or tags to be shown on the
#' plot.
#' @param weights A numeric vector giving the relative proportions of text size
#' corresponding to the given tag.
#' @param x,object An object of the type produced by tagcloud.
#' @param \dots Further arguments to be passed to downstream methods.
#' @param algorithm Name of the algorithm to use. Can be "oval", "fill",
#' "random", "snake", "list" or "clist". See Details.
#' @param scale If "auto", text expansion will be calculated automatically to
#' fit the available space. Otherwise, a numeric value used to modify the
#' calculated text sizes; tune scale to achieve a better fit.
#' @param scale.multiplier Multiplier for the final calculated text expansion
#' parameter. Increase if there is too much empty space around the tag cloud;
#' decrease if the tags go over the plot boundaries.
#' @param order Determines in which order the tags will be drawn. Can be
#' "size", "keep", "random", "height" or "width". See Details.
#' @param sel An integer or boolean vector indicating which terms from the
#' provided list will be used to plot. The vectors \code{col} and
#' \code{weights} will be filtered accordingly.
#' @param wmin All items in the \code{weights} vector smaller than \code{wmin}
#' will be changed to \code{wmin}
#' @param wmax All items in the \code{weights} vector larger than \code{wmax}
#' will be changed to \code{wmax}
#' @param floor Minimal text size. See Details.
#' @param ceiling Maximal text size. See Details.
#' @param family Font family to use, a vector containing font families to use
#' for each tag. For the \code{tagcloud} function, the special keyword "random"
#' can be used to assign random families (requires the extrafont package).
#' @param col Color or a vector containing colors to use for drawing the tags.
#' @param fvert Fraction of tags which will be rotated by 90 degrees
#' counterclockwise.
#' @param plot If FALSE, no plot will be produced.
#' @param add If TRUE, the tags will be added to the current plot instead of
#' creating a new plot.
#' @param with.box If TRUE, a rectangle will be plotted around each tag.
#' @return
#' 
#' \code{tagcloud} returns an object of the \code{tagcloud-class}, which really
#' is a data frame with the following columns:
#' 
#' \itemize{ 
#' \item\code{tags} -- the tags, words or phrases shown on the plot
#' \item\code{weights} -- a numeric vector that is used to calculate the size
#' of the plotted tags 
#' \item\code{family} -- name of the font family to be used in plotting 
#' \item\code{vertical} -- whether the tag should be rotated by 90
#' degrees counter-clockwise 
#' \item\code{x},\code{y} -- coordinates of the left
#' lower corner of the tags bounding box 
#' \item\code{w},\code{h} -- width and height of the bounding box 
#' \item\code{cex} -- text expansion factor, see \code{\link{par}} 
#' \item\code{s} -- surface of the tag (width x height) 
#' }
#' 
#' The object of the \code{tagcloud} class can be manipulated using
#' \code{\link{editor.tagcloud}} and displayed using \code{\link{plot}},
#' \code{\link{print}} and \code{\link{summary}} functions.
#' @note Care should be taken when using extra fonts loaded by the extrafont
#' package; not all fonts can be easily copied to a PDF file.
#' 
#' Some ideas in this package are based on the `wordcloud` package by Ian
#' Fellows.
#' @author January Weiner <january.weiner@@gmail.com>
#' @seealso \code{\link{editor.tagcloud}} -- interactive editing of tagcloud
#' objects.
#' 
#' \code{\link{strmultline}} -- splitting multi-word sentences into lines for a
#' better cloud display.
#' 
#' \code{\link{smoothPalette}} -- mapping values onto a color gradient.
#' @keywords tags word cloud weighted list tag cloud
#' @examples
#' 
#' 
#' # a rather boring tag cloud
#' data( gambia )
#' terms <- gambia$Term
#' tagcloud( terms )
#' 
#' # tag cloud with weights relative to P value
#' # colors relative to odds ratio, from light
#' # grey to black
#' weights <- -log( gambia$Pvalue )
#' colors  <- smoothPalette( gambia$OddsRatio, max=4 )
#' tagcloud( terms, weights, col= colors, algorithm= "oval" )
#' 
#' # tag cloud filling the whole plot
#' tagcloud( terms, weights, col= colors, algorithm= "fill" )
#' 
#' # just a list of only the first ten terms
#' tagcloud( terms, weights, sel= 1:10,
#'           col= colors, algorithm= "list", order= "width" )
#' 
#' # oval, with line breaks in terms
#' terms <- strmultline( gambia$Term )
#' tagcloud( terms, weights, col= colors, algorithm= "oval" )
#' 
#' \dontrun{
#' # shows available font families, scaled according to
#' # the total disk space occupied by the fonts
#' require( extrafont )
#' ft <- fonttable()
#' fi <- file.info( fonttable()$fontfile )
#' families <- unique( ft$FamilyName )
#' sizes    <- sapply( families,function( x ) sum( fi[ ft$FamilyName == x, "size" ] ) )
#' tagcloud( families, sizes, family= families )
#' }
#' 
#' 
#' @useDynLib tagcloud
#' @import Rcpp
#' @import RColorBrewer
#' @importFrom graphics par strwidth strheight identify locator plot plot.new plot.window rect text
#' @importFrom stats rnorm runif
#' @importFrom grDevices colorRampPalette
#' @export tagcloud
tagcloud <- function( tags, weights= 1, 
  algorithm= "oval", scale= "auto", scale.multiplier= 1,
  order= "size", sel= NULL,
  wmin= NULL, wmax= NULL, floor= 1, ceiling= 3, 
  family= NULL, col= NULL, 
  fvert= 0,
  plot= TRUE, add= FALSE
   ) {

  tags   <- as.character( tags )
  n.tags <- length( tags )

  # meta holds meta information about the tags: colors, weights and font
  meta <- data.frame( tags= tags )
  meta$weights <- weights
  meta$family  <- family
  meta$colors  <- col

  # boxes hold the actual positioning of the tags on the plot
  boxes <- matrix( 0, nrow= n.tags, ncol= 7 )
  colnames( boxes ) <- c( "x", "y", "w", "h", "cex", "s", "srt" )

  # random font family specification - needs extrafont
  if( ! missing( family ) && length( family ) == 1 && family == "random" ) {
    if(requireNamespace( "extrafont", quietly=TRUE )) {
      meta$family <- sample( extrafont::fonts(), n.tags )
    } else {
      warning("Install the 'extrafont' package for additional fonts")
    }
  }

  if( ! missing( sel ) ) {
    meta   <- meta[ sel, ] 
    boxes  <- boxes[ sel, ]
    n.tags <- nrow( boxes ) 
  }

  if( plot & ! add ) {
    plot.new()
    old.par <- par( mar= c( 0, 0, 0, 0 ) )
    plot.window( xlim= c( 0, 1 ), ylim= c( 0, 1 ), asp= 1 )
  }

  boxes[, "cex"] <- .calc.cex( meta$weights, floor, ceiling, wmin= wmin, wmax= wmax )
  boxes <- .calc.sizes( meta$tags, boxes, meta$family )

  if( scale == "auto" ) scale <- .auto.scale( boxes )
  scale <- scale * scale.multiplier

  boxes[,"cex"] <- boxes[,"cex"] * scale
  boxes <- .calc.sizes( meta$tags, boxes, meta$family )

  meta$vertical <- rep( FALSE, n.tags )
  if( length( fvert ) > 1 ) {
    meta$vertical <- fvert
  } else if( fvert > 0 ) {
    n <- as.integer( n.tags * fvert )
    #meta$vertical[ order( boxes[,"w"] )[ 1:n ] ] <- TRUE 
    meta$vertical[ runif( n.tags, 0, 1 ) < fvert ] <- TRUE
  }

  order <- match.arg( order, c( "size", "keep", "random", "width", "height" ) )
  order <- switch( order,
    keep= 1:nrow( boxes ),
    random= sample( 1:nrow( boxes ) ),
    size= order( boxes[,"s"], decreasing= T ),
    width= order( boxes[,"w"], decreasing= T ),
    height= order( boxes[,"h"], decreasing= T )
     )

  meta  <- meta[ order, ]
  boxes <- boxes[ order, ]

  # filter out invisible tags
  sel   <- boxes[,"h"] < 1e-6 | boxes[,"w"] < 1e-6
  meta  <- meta[ ! sel, ]
  boxes <- boxes[ ! sel, ]

  boxes[ meta$vertical, c( "w", "h" ) ] <- boxes[ meta$vertical, c( "h", "w" ) ]

  algorithm <- match.arg( algorithm, 
    c( "oval", "fill", "random", "list", "clist", "snake" ) )

  boxes <- switch( algorithm,
      snake=algorithm.snake( boxes, meta ),
      oval=algorithm.spiral( boxes ),
      normal=algorithm.normal( boxes ),
      fill=algorithm.ulam( boxes ),
      random=algorithm.random( boxes ),
      clist=algorithm.list( boxes, meta, centered= TRUE ), 
      list=algorithm.list( boxes, meta ) )

  debugprf( "final surface index: %.3f", boxes.ratio( boxes ) )

  boxes <- .tagcloud.new.object( boxes, meta, algorithm, scale )

  if( plot )  {
    debugprf( "plotting" )
    plot( boxes, add= T, with.box= tagcloud.debug )
    if( ! add ) par( old.par ) 
  }

  return( invisible( boxes ))
}
