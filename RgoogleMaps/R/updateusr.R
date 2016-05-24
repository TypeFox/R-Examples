`updateusr` <-structure(function#Updates the 'usr' coordinates in the current plot.
### For a traditional graphics plot this function will update the 'usr'
### coordinates by transforming a pair of points from the current usr
### coordinates to those specified.
(
  x1, ##<< The x-coords of 2 points in the current 'usr' coordinates, or anything that can be passed to \code{xy.coords}.
  y1 = NULL, ##<< The y-coords of 2 points in the current 'usr' coordinates, or an object representing the points in the new 'usr' coordinates. 
  x2, ##<< The x-coords for the 2 points in the new coordinates.
  y2 = NULL ##<< The y-coords for the 2 points in the new coordinates.
) 
{
    xy1 <- xy.coords(x1, y1)
    xy2 <- if (missing(x2) && missing(y2)) {
        xy.coords(y1)
    }
    else {
        xy.coords(x2, y2)
    }
    cur.usr <- par("usr")
    xslope <- diff(xy2$x)/diff(xy1$x)
    yslope <- diff(xy2$y)/diff(xy1$y)
    new.usr.x <- xslope * (cur.usr[1:2] - xy1$x) + xy2$x
    new.usr.y <- yslope * (cur.usr[3:4] - xy1$y) + xy2$y
    invisible(par(usr = c(new.usr.x, new.usr.y)))

##details<< Sometimes graphs (in the traditional graphing scheme) end up with usr
##  coordinates different from expected for adding to the plot (for
##  example \code{barplot} does not center the bars at integers).  This
##  function will take 2 points in the current 'usr' coordinates and the
##  desired 'usr' coordinates of the 2 points and transform the user
##  coordinates to make this happen.  The updating only shifts and scales
##  the coordinates, it does not do any rotation or warping transforms.
##
##  If \code{x1} and \code{y1} are lists or matricies and \code{x2} and
##  \code{y2} are not specified, then \code{x1} is taken to be the
##  coordinates in the current system and \code{y1} is the coordinates in
##  the new system.
##
##  Currently you need to give the function exactly 2 points in each
##  system.  The 2 points cannot have the same x values or y values in
##  either system.

##note<<  Currently you need to give coordinates for exactly 2 points without
##  missing values.  Future versions of the function will allow missing
##  values or multiple points.
##  
##  Note by Markus Loecher: both the source and the documentations were copied from the package TeachingDemos version 2.3


### An invisible list with the previous 'usr' coordinates from \code{par}.
}, ex = function(){
  tmp <- barplot(1:4)
  updateusr(tmp[1:2], 0:1, 1:2, 0:1)
  lines(1:4, c(1,3,2,2), lwd=3, type='b',col='red')

  # update the y-axis to put a reference distribution line in the bottom
  # quarter

  tmp <- rnorm(100)
  hist(tmp)
  tmp2 <- par('usr')
  xx <- seq(min(tmp), max(tmp), length.out=250)
  yy <- dnorm(xx, mean(tmp), sd(tmp))
  updateusr( tmp2[1:2], tmp2[3:4], tmp2[1:2], c(0, max(yy)*4) )
  lines(xx,yy)

})

