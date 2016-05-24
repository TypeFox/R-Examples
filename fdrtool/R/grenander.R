### grenander.R  (2007-06-13)
###
###     Grenander Density Estimator
###
### Copyright 2006-2007 Korbinian Strimmer 
###
###
### This file is part of the `fdrtool' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


grenander = function(F, type=c("decreasing", "increasing"))
{
  if( !any(class(F) == "ecdf") ) stop("ecdf object required as input!")

  type <- match.arg(type)

  if (type == "decreasing")
  {
    # find least concave majorant of ECDF
    ll = gcmlcm(environment(F)$x, environment(F)$y, type="lcm")
  }
  else
  {
    # find greatest convex minorant of ECDF
    l = length(environment(F)$y)
    ll = gcmlcm(environment(F)$x, c(0,environment(F)$y[-l]), type="gcm")
  }

  f.knots = ll$slope.knots
  f.knots = c(f.knots, f.knots[length(f.knots)])

  g = list(F=F,
       x.knots=ll$x.knots,
       F.knots=ll$y.knots,
       f.knots=f.knots)

  class(g) <- "grenander"
  
  return(g)
}

plot.grenander <- function(x, ...)
{
  if (x$f.knots[1] > x$f.knots[2])
    main = "Grenander Decreasing Density"
  else
    main = "Grenander Increasing Density"

  par(mfrow=c(1,2))

  plot(x$x.knots, x$f.knots, type="s", xlab="x", ylab="fn(x)",
     main=main, col=4, lwd=2, ...)
 
  plot(x$F, do.points=FALSE)
  lines(x$x.knots, x$F.knots, type='l', col=4, lwd=2)


  par(mfrow=c(1,1))
}



