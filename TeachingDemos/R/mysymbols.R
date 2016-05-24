my.symbols <- function(x, y=NULL, symb, inches=1, xsize, ysize,
                       add=TRUE,
                       vadj=0.5, hadj=0.5,
                       symb.plots=FALSE,
                       xlab=deparse(substitute(x)),
                       ylab=deparse(substitute(y)), main=NULL,
                       xlim=NULL, ylim=NULL, linesfun=lines,
                       ..., MoreArgs ) {

  if(!add){
	plot(x,y, type='n', xlab=xlab,ylab=ylab,
             xlim=xlim,ylim=ylim,main=main)
  }

  xy <- xy.coords(x,y,recycle=TRUE)

  pin <- par('pin')
  usr <- par('usr')
  usr.x <- usr[2] - usr[1]
  usr.y <- usr[4] - usr[3]

#  tmp <- cnvrt.coords(xy,input='usr')$plt
  tmp <- list()
  tmp$x <- grconvertX(xy$x, to='npc')
  tmp$y <- grconvertY(xy$y, to='npc')

  tmp.xlen <- length(tmp$x)

  if( (length(inches) != 1) && (length(inches) != tmp.xlen) ) {
    inches <- rep(inches, length.out=tmp.xlen)
  }
  if( (length(hadj) != 1) && (length(hadj) != tmp.xlen) ) {
    hadj <- rep(hadj, length.out=tmp.xlen)
  }
  if( (length(vadj) != 1) && (length(vadj) != tmp.xlen) ) {
    vadj <- rep(vadj, length.out=tmp.xlen)
  }

  if( missing(xsize) ) {
      if( missing(ysize) ) { # use inches
          x.low  <- tmp$x -    hadj *inches/pin[1]
          x.high <- tmp$x + (1-hadj)*inches/pin[1]
          y.low  <- tmp$y -    vadj *inches/pin[2]
          y.high <- tmp$y + (1-vadj)*inches/pin[2]
      } else { # ysize only
          y.low  <- tmp$y - vadj*ysize/usr.y
          y.high <- tmp$y + (1-vadj)*ysize/usr.y
          x.low  <- tmp$x - hadj/pin[1]*pin[2]/usr.y*ysize
          x.high <- tmp$x + (1-hadj)/pin[1]*pin[2]/usr.y*ysize
      }
  } else {
      if( missing(ysize) ) { # xsize only
          x.low  <- tmp$x - hadj*xsize/usr.x
          x.high <- tmp$x + (1-hadj)*xsize/usr.x
          y.low  <- tmp$y - vadj/pin[2]*pin[1]/usr.x*xsize
          y.high <- tmp$y + (1-vadj)/pin[2]*pin[1]/usr.x*xsize
      } else {  # both xsize and ysize
          x.low  <- tmp$x - hadj*xsize/usr.x
          x.high <- tmp$x + (1-hadj)*xsize/usr.x
          y.low  <- tmp$y - vadj*ysize/usr.y
          y.high <- tmp$y + (1-vadj)*ysize/usr.y
      }
  }


#  xy.low  <- cnvrt.coords(x.low,  y.low,  'plt')$fig
#  xy.high <- cnvrt.coords(x.high, y.high, 'plt')$fig

  xy.low <- list()
  xy.low$x <- grconvertX(x.low, from='npc', to='nfc')
  xy.low$y <- grconvertY(y.low, from='npc', to='nfc')

  xy.high <- list()
  xy.high$x <- grconvertX(x.high, from='npc', to='nfc')
  xy.high$y <- grconvertY(y.high, from='npc', to='nfc')


  plotfun <- if( is.function(symb) ) {
    if(symb.plots) {
      function(xlow,xhigh,ylow,yhigh,symb, ...) {
        op <- par(c('plt','usr','xpd'))
        on.exit(par(op))
        par(xpd=TRUE)
        par(plt=c(xlow,xhigh,ylow,yhigh), new=TRUE)
        par(usr=c(-1,1,-1,1))
        symb(...)
      }
    } else {
      function(xlow,xhigh,ylow,yhigh,symb, ...) {
        op <- par(c('plt','usr','xpd'))
        on.exit(par(op))
        par(xpd=TRUE)
        par(plt=c(xlow,xhigh,ylow,yhigh))
        par(usr=c(-1,1,-1,1))
        suppressWarnings(
            linesfun( symb(...), ... )
                       )
      }
    }
  } else {
    function(xlow,xhigh,ylow,yhigh,symb, ...) {
      op <- par(c('plt','usr','xpd'))
      on.exit(par(op))
      par(xpd=TRUE)
      par(plt=c(xlow,xhigh,ylow,yhigh))
      par(usr=c(-1,1,-1,1))
      linesfun(symb, ...)
    }
  }

  funargs <- list(xlow=xy.low$x, xhigh=xy.high$x,
                        ylow=xy.low$y, yhigh=xy.high$y)
  if( length(list(...)) ) {
    funargs <- c(funargs,
                 lapply(list(...), function(x) rep(x,length.out=tmp.xlen) )
                 )
  }

  funargs$FUN <- plotfun
  if (missing(MoreArgs)) {
    funargs$MoreArgs <- list(symb=symb)
  } else {
    funargs$MoreArgs <- c(MoreArgs, list(symb=symb))
  }

  do.call(mapply, funargs)

  invisible(NULL)
}



ms.male <- structure(c(0, 0.022, 0.0439, 0.0657, 0.0874, 0.109, 0.1303,
0.1514, 0.1722, 0.1926, 0.2127, 0.2324, 0.2516, 0.2703, 0.2885,
0.3062, 0.3233, 0.3397, 0.3555, 0.3706, 0.385, 0.3986, 0.4115,
0.4236, 0.4348, 0.4453, 0.4548, 0.4635, 0.4713, 0.4782, 0.4841,
0.4891, 0.4932, 0.4964, 0.4985, 0.4997, 0.5, 0.4992, 0.4976,
0.4949, 0.4913, 0.4868, 0.4813, 0.4748, 0.4675, 0.4593, 0.4501,
0.4401, 0.4293, 0.4176, 0.4052, 0.3919, 0.3779, 0.3631, 0.3477,
0.3316, 0.3148, 0.2974, 0.2795, 0.261, 0.242, 0.2226, 0.2027,
0.1824, 0.1618, 0.1409, 0.1197, 0.0982, 0.0766, 0.0548, 0.0329,
0.011, -0.011, -0.0329, -0.0548, -0.0766, -0.0982, -0.1197, -0.1409,
-0.1618, -0.1824, -0.2027, -0.2226, -0.242, -0.261, -0.2795,
-0.2974, -0.3148, -0.3316, -0.3477, -0.3631, -0.3779, -0.3919,
-0.4052, -0.4176, -0.4293, -0.4401, -0.4501, -0.4593, -0.4675,
-0.4748, -0.4813, -0.4868, -0.4913, -0.4949, -0.4976, -0.4992,
-0.5, -0.4997, -0.4985, -0.4964, -0.4932, -0.4891, -0.4841, -0.4782,
-0.4713, -0.4635, -0.4548, -0.4453, -0.4348, -0.4236, -0.4115,
-0.3986, -0.385, -0.3706, -0.3555, -0.3397, -0.3233, -0.3062,
-0.2885, -0.2703, -0.2516, -0.2324, -0.2127, -0.1926, -0.1722,
-0.1514, -0.1303, -0.109, -0.0874, -0.0657, -0.0439, -0.022,
0, NA, 0.3536, 1, 0.6, NA, 1, 1, 0.5, 0.4995, 0.4981, 0.4957,
0.4923, 0.488, 0.4827, 0.4765, 0.4694, 0.4614, 0.4525, 0.4427,
0.4321, 0.4206, 0.4083, 0.3953, 0.3814, 0.3669, 0.3516, 0.3357,
0.3191, 0.3018, 0.284, 0.2657, 0.2468, 0.2275, 0.2077, 0.1875,
0.167, 0.1461, 0.125, 0.1036, 0.082, 0.0603, 0.0384, 0.0165,
-0.0055, -0.0274, -0.0494, -0.0712, -0.0928, -0.1143, -0.1356,
-0.1566, -0.1773, -0.1977, -0.2176, -0.2372, -0.2563, -0.2749,
-0.293, -0.3105, -0.3274, -0.3437, -0.3593, -0.3743, -0.3885,
-0.4019, -0.4146, -0.4265, -0.4375, -0.4477, -0.4571, -0.4655,
-0.4731, -0.4797, -0.4855, -0.4903, -0.4941, -0.497, -0.4989,
-0.4999, -0.4999, -0.4989, -0.497, -0.4941, -0.4903, -0.4855,
-0.4797, -0.4731, -0.4655, -0.4571, -0.4477, -0.4375, -0.4265,
-0.4146, -0.4019, -0.3885, -0.3743, -0.3593, -0.3437, -0.3274,
-0.3105, -0.293, -0.2749, -0.2563, -0.2372, -0.2176, -0.1977,
-0.1773, -0.1566, -0.1356, -0.1143, -0.0928, -0.0712, -0.0494,
-0.0274, -0.0055, 0.0165, 0.0384, 0.0603, 0.082, 0.1036, 0.125,
0.1461, 0.167, 0.1875, 0.2077, 0.2275, 0.2468, 0.2657, 0.284,
0.3018, 0.3191, 0.3357, 0.3516, 0.3669, 0.3814, 0.3953, 0.4083,
0.4206, 0.4321, 0.4427, 0.4525, 0.4614, 0.4694, 0.4765, 0.4827,
0.488, 0.4923, 0.4957, 0.4981, 0.4995, 0.5, NA, 0.3536, 1, 1,
NA, 1, 0.6), .Dim = as.integer(c(151, 2)))

ms.female <- structure(c(0, 0.022, 0.0439, 0.0657, 0.0874, 0.109, 0.1303,
0.1514, 0.1722, 0.1926, 0.2127, 0.2324, 0.2516, 0.2703, 0.2885,
0.3062, 0.3233, 0.3397, 0.3555, 0.3706, 0.385, 0.3986, 0.4115,
0.4236, 0.4348, 0.4453, 0.4548, 0.4635, 0.4713, 0.4782, 0.4841,
0.4891, 0.4932, 0.4964, 0.4985, 0.4997, 0.5, 0.4992, 0.4976,
0.4949, 0.4913, 0.4868, 0.4813, 0.4748, 0.4675, 0.4593, 0.4501,
0.4401, 0.4293, 0.4176, 0.4052, 0.3919, 0.3779, 0.3631, 0.3477,
0.3316, 0.3148, 0.2974, 0.2795, 0.261, 0.242, 0.2226, 0.2027,
0.1824, 0.1618, 0.1409, 0.1197, 0.0982, 0.0766, 0.0548, 0.0329,
0.011, -0.011, -0.0329, -0.0548, -0.0766, -0.0982, -0.1197, -0.1409,
-0.1618, -0.1824, -0.2027, -0.2226, -0.242, -0.261, -0.2795,
-0.2974, -0.3148, -0.3316, -0.3477, -0.3631, -0.3779, -0.3919,
-0.4052, -0.4176, -0.4293, -0.4401, -0.4501, -0.4593, -0.4675,
-0.4748, -0.4813, -0.4868, -0.4913, -0.4949, -0.4976, -0.4992,
-0.5, -0.4997, -0.4985, -0.4964, -0.4932, -0.4891, -0.4841, -0.4782,
-0.4713, -0.4635, -0.4548, -0.4453, -0.4348, -0.4236, -0.4115,
-0.3986, -0.385, -0.3706, -0.3555, -0.3397, -0.3233, -0.3062,
-0.2885, -0.2703, -0.2516, -0.2324, -0.2127, -0.1926, -0.1722,
-0.1514, -0.1303, -0.109, -0.0874, -0.0657, -0.0439, -0.022,
0, NA, 0, 0, NA, -0.25, 0.25, 0.5, 0.4995, 0.4981, 0.4957, 0.4923,
0.488, 0.4827, 0.4765, 0.4694, 0.4614, 0.4525, 0.4427, 0.4321,
0.4206, 0.4083, 0.3953, 0.3814, 0.3669, 0.3516, 0.3357, 0.3191,
0.3018, 0.284, 0.2657, 0.2468, 0.2275, 0.2077, 0.1875, 0.167,
0.1461, 0.125, 0.1036, 0.082, 0.0603, 0.0384, 0.0165, -0.0055,
-0.0274, -0.0494, -0.0712, -0.0928, -0.1143, -0.1356, -0.1566,
-0.1773, -0.1977, -0.2176, -0.2372, -0.2563, -0.2749, -0.293,
-0.3105, -0.3274, -0.3437, -0.3593, -0.3743, -0.3885, -0.4019,
-0.4146, -0.4265, -0.4375, -0.4477, -0.4571, -0.4655, -0.4731,
-0.4797, -0.4855, -0.4903, -0.4941, -0.497, -0.4989, -0.4999,
-0.4999, -0.4989, -0.497, -0.4941, -0.4903, -0.4855, -0.4797,
-0.4731, -0.4655, -0.4571, -0.4477, -0.4375, -0.4265, -0.4146,
-0.4019, -0.3885, -0.3743, -0.3593, -0.3437, -0.3274, -0.3105,
-0.293, -0.2749, -0.2563, -0.2372, -0.2176, -0.1977, -0.1773,
-0.1566, -0.1356, -0.1143, -0.0928, -0.0712, -0.0494, -0.0274,
-0.0055, 0.0165, 0.0384, 0.0603, 0.082, 0.1036, 0.125, 0.1461,
0.167, 0.1875, 0.2077, 0.2275, 0.2468, 0.2657, 0.284, 0.3018,
0.3191, 0.3357, 0.3516, 0.3669, 0.3814, 0.3953, 0.4083, 0.4206,
0.4321, 0.4427, 0.4525, 0.4614, 0.4694, 0.4765, 0.4827, 0.488,
0.4923, 0.4957, 0.4981, 0.4995, 0.5, NA, -0.5, -1, NA, -0.8,
-0.8), .Dim = as.integer(c(150, 2)))

ms.polygon <- function(n, r=1, adj=pi/2, ...) {
  tmp <- seq(0,2*pi, length.out=n+1) + adj
  cbind(cos(tmp), sin(tmp)) * r
}

ms.filled.polygon <- function(n, r=1, adj=pi/2, fg=par('fg'),
                              bg=par('fg'), ... ) {
  tmp <- seq(0,2*pi, length.out=n+1) + adj
  polygon(cos(tmp)*r,sin(tmp)*r, col=bg, border=fg, ...)
  NULL
}


ms.polygram <- function(n, r=1, adj=pi/2, ...) {
  if (n == 1) {
    return(rbind( c(0,0), c(cos(adj),sin(adj))*r))
  }
  if (n == 2) {
    return(rbind( c(cos(adj),sin(adj)), c(cos(adj+pi),sin(adj+pi))) * r)
  }
  if (n == 3) {
    return(rbind(
                 c(0,0),
                 c(cos(adj),sin(adj)),
                 NA,
                 c(0,0),
                 c(cos(adj+2*pi/3), sin(adj+2*pi/3)),
                 NA,
                 c(0,0),
                 c(cos(adj+4*pi/3), sin(adj+4*pi/3)))*r)
  }
  if (n == 4) {
    return(rbind(
                 c(cos(adj),sin(adj)),
                 c(cos(adj+pi),sin(adj+pi)),
                 NA,
                 c(cos(adj+pi/2), sin(adj+pi/2)),
                 c(cos(adj+3*pi/2), sin(adj+3*pi/2))) * r )
  }
  if (n == 6) {
    tmp <- c( 0, 2*pi/3, 4*pi/3, 2*pi )
    tmp <- c(tmp, NA, tmp+pi/3)+adj
    return( cbind( cos(tmp), sin(tmp) )*r )
  }
  skp <- floor( n/2 - 0.1 )
  tmp <- seq( 0, skp*2*pi, length.out=n+1 ) + adj
  tmp2 <- cbind(cos(tmp), sin(tmp))*r
  while( any( duplicated( round( tmp2[-1,], 5 ) ) ) ){
    skp <- skp - 1
    tmp <- seq( 0, skp*2*pi, length.out=n+1 ) + adj
    tmp2 <- cbind( cos(tmp), sin(tmp))*r
  }
  return(tmp2)
}

ms.arrows <- function(angle, r=1, adj=0.5, length=0.1, ...) {
  fr <- c( cos(angle), sin(angle) ) * (-r) * adj
  to <- c( cos(angle), sin(angle) ) * r * (1-adj)
  arrows(fr[1],fr[2],to[1],to[2], length=length, ...)
  NULL
}

ms.sunflowers <- function(n,r=0.3,adj=pi/2, ...) {
  tmp <- seq(0,2*pi, length.out=36)
  tmp2 <- cbind( cos(tmp), sin(tmp) ) * r
  tmp <- seq( 0, 2*pi, length.out=n+1 )[-(n+1)] + adj
  tmp.x <- c(rbind(NA,cos(tmp)*r, cos(tmp)))
  tmp.y <- c(rbind(NA,sin(tmp)*r, sin(tmp)))
  rbind(tmp2, cbind(tmp.x, tmp.y) )
}


ms.image <- function(img, transpose=TRUE, ...) {
    d <- dim(img)

    cols <- if(d[3] == 3) {
        rgb(img[,,1], img[,,2], img[,,3])
    } else if(d[3] == 4) {
        rgb(img[,,1], img[,,2], img[,,3], img[,,4])
    } else {
        stop('image must be array with 3rd dimension equal to 3 or 4')
    }

    if(transpose) {
        tmp <- matrix( seq(length=d[1]*d[2]), ncol=d[1], byrow=TRUE)
        tmp <- tmp[ , rev(seq(length=d[1])) ]
    } else {
        tmp <- matrix( seq(length=d[1]*d[2]), ncol=d[2] )
        tmp <- tmp[ , rev(seq(length=d[2])) ]
    }

    image(tmp, col=cols, axes=FALSE, xlab='', ylab='')
}



# do ms.raster that uses rasterImage





