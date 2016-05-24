tbox <-
function( txt, x, y, wd, ht,
          font=2, lwd=2,
          col.txt=par("fg"),
          col.border=par("fg"),
          col.bg="transparent" )
{
rect( x-wd/2, y-ht/2, x+wd/2, y+ht/2,
      lwd=lwd, border=col.border, col=col.bg )
text( x, y, txt, font=font, col=col.txt )
invisible( c( x, y, wd, ht ) )
}

dbox <-
function( x, y, wd, ht=wd,
          font=2, lwd=2, cwd=5,
          col.cross=par("fg"),
          col.border=par("fg"),
          col.bg="transparent" )
{
rect( x-wd/2, y-ht/2, x+wd/2, y+ht/2, lwd=lwd, border=col.border, col=col.bg )
ch <- ht*2/3
segments( c(x     , x-ch/3),
          c(y+ch/2, y+ch/6),
          c(x     , x+ch/3),
          c(y-ch/2, y+ch/6), lwd=cwd, col=col.cross )
invisible( c( x, y, wd, ht ) )
}

fillarr <-
function( x1, y1, x2, y2, gap=2, fr=0.8,
          angle=17, lwd=2, length=par("pin")[1]/30, ... )
{
fr <- 1-gap/sqrt((x1-x2)^2+(y1-y2)^2)
if( !missing(fr)  ) if( fr > 1 ) fr <- fr/100
for( a in 1:angle )
arrows( x1 + (x2-x1)*(1-fr)/2,
        y1 + (y2-y1)*(1-fr)/2,
        x2 - (x2-x1)*(1-fr)/2,
        y2 - (y2-y1)*(1-fr)/2, angle=a, lwd=lwd, ... )
}

std.vec <-
function( a, b )
{
  l <- sqrt(a^2+b^2)
  if( l==0 )
  return( c(0,0) )
  else
  return( c(a/l,b/l) )
}

boxarr <-
function (b1, b2, offset = FALSE, pos = 0.45, ...)
{
    d <- std.vec(b2[1] - b1[1], b2[2] - b1[2])
    dd <- d * offset
    x1 <- b1[1] - dd[2]
    y1 <- b1[2] + dd[1]
    w1 <- b1[3]
    h1 <- b1[4]
    x2 <- b2[1] - dd[2]
    y2 <- b2[2] + dd[1]
    w2 <- b2[3]
    h2 <- b2[4]
    hx1 <- x1 + ifelse((y2-y1) != 0, (x2-x1) * ((h1/2)/abs(y2-y1)), sign(x2-x1) * w1/2)
    vx1 <- x1 + ifelse((x2-x1) != 0, (x2-x1) * ((w1/2)/abs(x2-x1)), 0)
    hx2 <- x2 + ifelse((y1-y2) != 0, (x1-x2) * ((h2/2)/abs(y1-y2)), sign(x1-x2) * w2/2)
    vx2 <- x2 + ifelse((x1-x2) != 0, (x1-x2) * ((w2/2)/abs(x1-x2)), 0)
    hy1 <- y1 + ifelse((y2-y1) != 0, (y2-y1) * ((h1/2)/abs(y2-y1)), 0)
    vy1 <- y1 + ifelse((x2-x1) != 0, (y2-y1) * ((w1/2)/abs(x2-x1)), sign(y2-y1) * h1/2)
    hy2 <- y2 + ifelse((y1-y2) != 0, (y1-y2) * ((h2/2)/abs(y1-y2)), 0)
    vy2 <- y2 + ifelse((x1-x2) != 0, (y1-y2) * ((w2/2)/abs(x1-x2)), sign(y1-y2) * h2/2)
if( abs(vy1-y1) < h1/2 ) { bx1 <- vx1
                           by1 <- vy1 }
                    else { bx1 <- hx1
                           by1 <- hy1 }
if( abs(vy2-y2) < h2/2 ) { bx2 <- vx2
                           by2 <- vy2 }
                    else { bx2 <- hx2
                           by2 <- hy2 }
fillarr( bx1, by1, bx2, by2, ... )
invisible( list( x = bx1*(1-pos)+bx2*pos,
                 y = by1*(1-pos)+by2*pos,
                 d = d ) )
}

wh.no <-
function( tt, i, j )
{
## Utility to count the number of non-NA off diagonal elements with
## row<i or ( col=i & col<=j )
diag(tt) <- NA
sum(!is.na(tt[1:i,])) - sum(!is.na(tt[i,])) + sum(!is.na(tt[i,1:j]))
}

boxes <- function (obj, ...) UseMethod("boxes")

boxes.matrix <-
function( obj, ... )
  {
  Epi::boxes.Lexis( obj, ... )
  }

boxes.Lexis <-
function( obj, boxpos = FALSE,
                wmult = 1.15,
                hmult = 1.15,
                  cex = 1.45,
               show   = inherits( obj, "Lexis" ),
               show.Y = show,
              scale.Y = 1,
             digits.Y = 1,
              show.BE = FALSE,
               BE.sep = c("","","          ",""),
               show.D = show,
              scale.D = FALSE,
             digits.D = as.numeric(as.logical(scale.D)),
               show.R = is.numeric(scale.R),
              scale.R = 1,
             digits.R = as.numeric(as.logical(scale.R)),
               DR.sep = if( show.D ) c("\n(",")") else c("",""),
                eq.wd = TRUE,
                eq.ht = TRUE,
                   wd,
                   ht,
               subset = NULL,
              exclude = NULL,
                 font = 2,
                  lwd = 2,
           col.txt    = par("fg"),
           col.border = col.txt,
           col.bg     = "transparent",
           col.arr    = par("fg"),
              lwd.arr = 2,
             font.arr = 2,
              pos.arr = 0.45,
              txt.arr = NULL,
          col.txt.arr = col.arr,
           offset.arr = 2, ... )
{
if( inherits(obj,"Lexis") )
  {
  if( !is.factor(obj$lex.Cst) | !is.factor(obj$lex.Xst) ) obj <- factorize( obj )
  tm <- tmat( obj, Y=TRUE )
  tt <- tmat( obj, Y=FALSE )
  # If show.BE is a character, check if zeros should be printed
  noz <- FALSE
  if( is.character( show.BE ) )
    {
    noz <- ( tolower(show.BE) %in% c("nz","noz","nozero","nozeros") )
    show.BE <- TRUE
    }
  if( show.BE )
    {
    # Derive the persons at start and at end of study in each state
    Beg <- table( status( obj, at="entry", by.id=TRUE ) )
    End <- table( status( obj, at="exit" , by.id=TRUE ) )
    BE.line <- paste( ifelse( noz & Beg==0,
                              rep("  ",nchar(BE.sep[1])+2),
                              paste(BE.sep[1],
                                    formatC( Beg,
                                             format="f",
                                             digits=0,
                                             big.mark="," ),
                                    BE.sep[2],sep="") ),
                      ifelse( noz & End==0,
                              rep("  ",nchar(BE.sep[2])+nchar(BE.sep[3])+2),
                              paste(BE.sep[3],
                                    formatC( End,
                                             format="f",
                                             digits=0,
                                             big.mark="," ),
                                    BE.sep[4],sep="") ),
                      sep="" )
    }
  }
else if( is.matrix(obj) & diff(dim(obj))==0 )
  {
  tm <- tt <- obj
  }
else stop( "First argument must be a Lexis object or a square matrix.\n" )

# Put the transitions into D and the diagonal into Y.
D <- tm
diag( D ) <- NA
Y <- diag( tm ) / scale.Y

# Explicitly given numbers to be put in boxes ?
if( is.numeric(show.Y) )
  {
  Y <- show.Y
  show.Y <- TRUE
  }

# Compute the rates - vectors are automatically expanded to matrices columnwise
R <- D / Y * ifelse(scale.R,scale.R,1)

# If no person-years available anywhere, they or rates cannot be shown
if( all(is.na(Y)) ) show.Y <- show.R <- FALSE

# Derive state names, no. states and no. transitions
                      st.nam <- colnames( tm )
if( is.null(st.nam) ) st.nam <- paste(1:ncol(tm))
            pl.nam <- st.nam
      n.st <- length( st.nam )
      n.tr <- sum( !is.na(tm) ) - sum( !is.na(diag(tm)) )

# No extra line with person-years when they are NA, put in no at
# beginnning / end of study
if( show.Y ) pl.nam <- gsub( "\\\nNA", "",
                             paste( st.nam,
                                    formatC( Y,
                                             format="f",
                                             digits=digits.Y,
                                             big.mark="," ),
                                    sep="\n" ) )
if( show.BE ) pl.nam <- paste( pl.nam, BE.line, sep="\n" )

# Any subsetting:
sbst <- 1:nrow(tm)
if( !is.null(exclude) )
  {
  if( is.character(exclude) )
    exclude <- match( exclude, rownames(tm) )
  sbst <- sbst[-exclude]
  }
if( !is.null(subset) )
  {
  if( is.character(subset) )
      sbst <- match( subset, rownames(tm) )
  else sbst <- subset
  }
subset <- sbst

# Recycling of box-arguments
if( !missing(ht) )
if( length(ht         )<n.st ) ht         <- rep(ht        ,n.st)[1:n.st]
if( !missing(wd) )
if( length(wd         )<n.st ) wd         <- rep(wd        ,n.st)[1:n.st]
if( length(font       )<n.st ) font       <- rep(font      ,n.st)[1:n.st]
if( length(lwd        )<n.st ) lwd        <- rep(lwd       ,n.st)[1:n.st]
if( length(col.border )<n.st ) col.border <- rep(col.border,n.st)[1:n.st]
if( length(col.bg     )<n.st ) col.bg     <- rep(col.bg    ,n.st)[1:n.st]
if( length(col.txt    )<n.st ) col.txt    <- rep(col.txt   ,n.st)[1:n.st]
# Recycling of arrow-arguments
if( length(col.arr    )<n.tr ) col.arr    <- rep(col.arr    ,n.tr)[1:n.tr]
if( length(col.txt.arr)<n.tr ) col.txt.arr<- rep(col.txt.arr,n.tr)[1:n.tr]
if( length(lwd.arr    )<n.tr ) lwd.arr    <- rep(lwd.arr    ,n.tr)[1:n.tr]
if( length(font.arr   )<n.tr ) font.arr   <- rep(font.arr   ,n.tr)[1:n.tr]
if( length(pos.arr    )<n.tr ) pos.arr    <- rep(pos.arr    ,n.tr)[1:n.tr]

# Here comes the plot
# First setting up the plot area, and restoring the plot parameters later
opar <- par( mar=c(0,0,0,0), cex=cex )
on.exit( par(opar) )
plot( NA, bty="n",
      xlim=0:1*100, ylim=0:1*100, xaxt="n", yaxt="n", xlab="", ylab="" )
# String height and width only meaningful after a plot has been called
if( missing(ht) )
  {
  ht <- strheight( pl.nam ) * hmult
  if( eq.ht ) ht <- rep( max(ht), length(ht) )
  }
if( missing(wd) )
  {
  wd <- strwidth(  pl.nam ) * wmult
  if( eq.wd ) wd <- rep( max(wd), length(wd) )
  }

# If not supplied, ask for positions of boxes
if( is.list(boxpos) )
  {
  names(boxpos) <- tolower( names(boxpos) )
  if( length(intersect(names(boxpos),c("x","y")))<2 )
    stop( "The list given in 'boxpos=' must have components 'x' and 'y'" )
  if( length(boxpos$x) != n.st | length(boxpos$y) != n.st )
    stop( "The elements 'x' and 'y' of boxpos must both have length equal to no. states", n.st )
  xx <- boxpos$x
  yy <- boxpos$y
  }
if( is.logical(boxpos) )
  if( boxpos )
  {
  ang <- pi - 2*pi*((1:n.st-0.5)/n.st)
  xx <- cos( ang ) * 35 + 50
  yy <- sin( ang ) * 35 + 50
  }
else
  {
  xx <- yy <- numeric(n.st)
  for( i in subset )
     {
     cat( "\nClick for level ", st.nam[i] )
     flush.console()
     pt <- locator(1)
     xx[i] <- pt$x
     yy[i] <- pt$y
     tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i], ... )
     }
  cat( "\n" )
  }

# Plot the boxes and record position and size
b <- list()
for( i in subset ) b[[i]] <- tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i],
                                font=font[i],
                                  lwd=lwd[i],
                          col.txt=col.txt[i],
                    col.border=col.border[i],
                            col.bg=col.bg[i] )

# Arrows and text on them
arrowtext <- character(0)
for( i in subset ) for( j in subset )
  {
  if( !is.na(tt[i,j]) & i!=j )
    {
    # Which number of arrow is currently processed?
    a <- wh.no( tt, i, j )
    arr <- boxarr( b[[i]], b[[j]],
                   offset=(!is.na(tt[j,i]))*offset.arr,
                   lwd=lwd.arr[a], col=col.arr[a], pos=pos.arr[a], ... )
    if( show.D | show.R )
      {
    if( show.D & D[i,j]>0 )
      arrowtext[a] <- formatC( D[i,j], format="f", digits=0, big.mark="," )
    else arrowtext[a] <- ""
    if( show.R & R[i,j]>0 )
      arrowtext[a] <- paste( if( !is.null(arrowtext[a]) )
                                   paste( arrowtext[a], DR.sep[1], sep="" ),
                           formatC( R[i,j], format="f", digits=digits.R, big.mark="," ),
                          if( length(DR.sep) > 1 ) DR.sep[2], sep="" )
      }
    else
    if( !is.null(txt.arr) )
      arrowtext[a] <- txt.arr[a]
    if( !is.null(arrowtext[a]) )
    text( arr$x-arr$d[2], arr$y+arr$d[1],
          arrowtext[a],
          adj=as.numeric(c(arr$d[2]>0,arr$d[1]<0)),
          font=font.arr[a], col=col.txt.arr[a] )
    }
  }
# Redraw the boxes with white background to remove any arrows
for( i in subset ) tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i],
                         lwd=lwd[i], col.bg=par("bg") )
# Then redraw the boxes again
for( i in subset ) tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i],
                                font=font[i],
                                  lwd=lwd[i],
                          col.txt=col.txt[i],
                    col.border=col.border[i],
                            col.bg=col.bg[i] )

# Finally create an object with all information to re-draw the display
MSboxes <- list( Boxes = data.frame( xx = xx,
                                     yy = yy,
                                     wd = wd,
                                     ht = ht,
                                   font = font,
                                    lwd = lwd,
                                col.txt = col.txt,
                             col.border = col.border,
                                 col.bg = col.bg, stringsAsFactors=FALSE ),
           State.names = pl.nam,
                  Tmat = tt,
                Arrows = data.frame( lwd.arr = lwd.arr,
                                     col.arr = col.arr,
                                     pos.arr = pos.arr,
                                 col.txt.arr = col.txt.arr,
                                    font.arr = font.arr,
                                  offset.arr = offset.arr, stringsAsFactors=FALSE ),
             Arrowtext = arrowtext )
class( MSboxes ) <- "MS"
invisible( MSboxes )
}

boxes.MS <-
function( obj, sub.st, sub.tr, cex=1.5, ... )
{
if( !inherits(obj,"MS") )
  stop( "You must supply an object of class 'MSboxes'" )

n.st <- nrow( obj$Boxes )
n.tr <- nrow( obj$Arrows )
if( missing(sub.st) ) sub.st <- 1:n.st
if( missing(sub.tr) ) sub.tr <- 1:n.tr

# First setting up the plot area, and restoring the plot parameters later
opar <- par( mar=c(0,0,0,0), cex=cex )
on.exit( par(opar) )
plot( NA, bty="n",
      xlim=0:1*100, ylim=0:1*100, xaxt="n", yaxt="n", xlab="", ylab="" )

# Exercise the subsets by putting the relevant colors to "transparent"
obj$Boxes[-sub.st,c("col.txt",
                    "col.border",
                    "col.bg")] <- "transparent"
obj$Arrows[-sub.tr,c("col.arr",
                     "col.txt.arr")] <- "transparent"
# Then draw the boxes
b <- list()
for( i in 1:n.st ) b[[i]] <- with( obj$Boxes,
                    tbox( obj$State.names[i], xx[i], yy[i], wd[i], ht[i],
                                font=font[i],
                                  lwd=lwd[i],
                          col.txt=col.txt[i],
                    col.border=col.border[i],
                            col.bg=col.bg[i] ) )

# and the arrows
for( i in 1:n.st ) for( j in 1:n.st )
   {
  if( !is.na(obj$Tmat[i,j]) & i!=j )
    {
    a <- wh.no( obj$Tmat, i, j )
    arr <- with( obj$Arrows,
           boxarr( b[[i]], b[[j]],
                   offset=(!is.na(obj$Tmat[j,i]))*offset.arr,
                   lwd=lwd.arr[a],
                   col=col.arr[a],
                   pos=pos.arr[a], ... ) )
    with( obj$Arrows,
    text( arr$x-arr$d[2], arr$y+arr$d[1], obj$Arrowtext[a],
          adj=as.numeric(c(arr$d[2]>0,arr$d[1]<0)),
         font=font.arr[a],
          col=col.txt.arr[a] ) )
    }
   }
# Redraw the boxes with "bg" background to remove any arrows crossing
for( i in sub.st ) with( obj$Boxes,
           tbox( obj$State.names[i], xx[i], yy[i], wd[i], ht[i],
                       font=font[i],
                         lwd=lwd[i],
                     col.txt=par("bg"),
                  col.border=par("bg"),
                      col.bg=par("bg") ) )
# Then redraw the boxes again
for( i in sub.st ) with( obj$Boxes,
                   tbox( obj$State.names[i], xx[i], yy[i], wd[i], ht[i],
                               font=font[i],
                                 lwd=lwd[i],
                         col.txt=col.txt[i],
                   col.border=col.border[i],
                           col.bg=col.bg[i] ) )
# Done!
invisible( NULL )
}
