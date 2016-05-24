## Copyright (C) 2000 Paul Kienzle
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## usage: yi = interp1(x, y, xi [, 'method' [, 'extrap']])
##
## Interpolate the f <- function( x) at the points xi. The sample end  y
## points x must be strictly monotonic.  If y is a matrix with
## length(x) rows, yi will be a matrix of size rows(xi) by columns(y),
## or its transpose if xi is a row vector.
##
## Method is one of:
## 'nearest': return nearest neighbour.
## 'linear': linear interpolation from nearest neighbours
## 'pchip': piece-wise cubic hermite interpolating polynomial
## 'cubic': cubic interpolation from four nearest neighbours
## 'spline': cubic spline interpolation--smooth first and second
##           derivatives throughout the curve
## ['*' method]: same as method, but assumes x is uniformly spaced
##               only uses x(1) and x(2); usually faster, never slower
##
## Method defaults to 'linear'.
##
## If extrap is the string 'extrap', then extrapolate values beyond
## the endpoints.  If extrap is a number, replace values beyond the
## endpoints with that number.  If extrap is missing, assume NaN.
##
## Example:
##    xf=[0:0.05:10]; yf = sin(2*pi*xf/5)
##    xp=[0:10];      yp = sin(2*pi*xp/5)
##    lin=interp1(xp,yp,xf)
##    spl=interp1(xp,yp,xf,'spline')
##    cub=interp1(xp,yp,xf,'cubic')
##    near=interp1(xp,yp,xf,'nearest')
##    plot(xf,yf,';original;',xf,lin,';linear;',xf,spl,';spline;',...
##         xf,cub,';cubic;',xf,near,';nearest;',xp,yp,'*;;')
##
## See also: interp

## 2000-03-25 Paul Kienzle
##    added 'nearest' as suggested by Kai Habel
## 2000-07-17 Paul Kienzle
##    added '*' methods and matrix y
##    check for proper table lengths
## 2002-01-23 Paul Kienzle
##    fixed extrapolation

##interp1 <- function( x, y, xi, method, extrap)  {
##
##  if nargin<3 || nargin>5
##    usage("yi = interp1(x, y, xi [, 'method' [, 'extrap']])")
##  } #if
##
##  if nargin < 4, 
##    method = 'linear'
##  else
##    method = tolower(method)
##  } #if
##
##  if nargin < 5
##    extrap = NaN
##  } #if
##
##  ## reshape matrices for convenience
##  x = x(:)
##  if size(y,1)==1, y=y(:); } #if
##  transposed = (size(xi,1)==1)
##  xi = xi(:)
##
##  ## determine sizes
##  nx = size(x,1)
##  [ny, nc] = size(y)
##  if (nx < 2 || ny < 2)
##     error ("interp1: table too short")
##  } #if
##
##  ## determine which values are out of range and set them to extrap,
##  ## unless extrap=='extrap' in which case, extrapolate them like we
##  ## should be doing in the first place.
##  minx = x(1)
##  if (method(1) == '*')
##     dx = x(2) - x(1)
##     maxx = minx + (ny-1)*dx
##  else
##     maxx = x(nx)
##  } #if
##  if ischar(extrap) && strcmp(extrap,"extrap")
##    range=1:size(xi,1)
##    yi = matrix(0, size(xi,1), size(y,2))
##  else
##    range = find(xi >= minx & xi <= maxx)
##    yi = extrap*matrix(1, size(xi,1), size(y,2))
##    if isempty(range), 
##      if transposed, yi = yi.'; } #if
##      return
##    } #if
##    xi = xi(range)
##  } #if
##
##  if strcmp(method, 'nearest')
##    idx = lookup(0.5*(x(1:nx-1)+x(2:nx)), xi)+1
##    yi(range,:) = y(idx,:)
##
##  elseif strcmp(method, '*nearest')
##    idx = floor((xi-minx)/dx+1.5)
##    yi(range,:) = y(idx,:)
##
##  elseif strcmp(method, 'linear')
##    ## find the interval containing the test point
##    idx = lookup(x[2:(nx-1)], xi) + 1
##              # 2:(n-1) so that anything beyond the ends
##              # gets dumped into an interval
##    ## use the endpoints of the interval to define a line
##    dy = y(2:ny,:) - y(1:ny-1,:)
##    dx = x(2:nx) - x(1:nx-1)
##    s = (xi - x(idx))./dx(idx)
##    yi(range,:) = s(:,matrix(1, 1,nc)).*dy(idx,:) + y(idx,:)
##
##  elseif strcmp(method, '*linear')
##    ## find the interval containing the test point
##    t = (xi - minx)/dx + 1
##    idx = floor(t)
##
##    ## use the endpoints of the interval to define a line
##    dy = [y(2:ny,:) - y(1:ny-1,:); matrix(0, 1,nc)]
##    s = t - idx
##    yi(range,:) = s(:,matrix(1, 1,nc)).*dy(idx,:) + y(idx,:)
##
##  elseif strcmp(method, 'pchip') || strcmp(method, '*pchip')
##    if (nx == 2) x = linspace(minx, maxx, ny); } #if
##    yi(range,:) = pchip(x, y, xi)
##
##  elseif strcmp(method, 'cubic')
##    if (nx < 4 || ny < 4)
##      error ("interp1: table too short")
##    } #if
##    idx = lookup(x(3:nx-2), xi) + 1
##
##    ## Construct cubic equations for each interval using divided
##    ## differences (computation of c and d don't use divided differences
##    ## but instead solve 2 equations for 2 unknowns). Perhaps
##    ## reformulating this as a lagrange polynomial would be more efficient.
##    i=1:nx-3
##    J = matrix(1, 1,nc)
##    dx = diff(x)
##    dx2 = x(i+1).^2 - x(i).^2
##    dx3 = x(i+1).^3 - x(i).^3
##    a=diff(y,3)./dx(i,J).^3/6
##    b=(diff(y(1:nx-1,:),2)./dx(i,J).^2 - 6*a.*x(i+1,J))/2
##    c=(diff(y(1:nx-2,:),1) - a.*dx3(:,J) - b.*dx2(:,J))./dx(i,J)
##    d=y(i,:) - ((a.*x(i,J) + b).*x(i,J) + c).*x(i,J)
##    yi(range,:) = ((a(idx,:).*xi(:,J) + b(idx,:)).*xi(:,J) ...
##         + c(idx,:)).*xi(:,J) + d(idx,:)
##
##  elseif strcmp(method, '*cubic')
##    if (nx < 4 || ny < 4)
##      error ("interp1: table too short")
##    } #if
##
##    ## From: Miloje Makivic 
##    ## http://www.npac.syr.edu/projects/nasa/MILOJE/final/node36.html
##    t = (xi - minx)/dx + 1
##    idx = max(min(floor(t), ny-2), 2)
##    t = t - idx
##    t2 = t.*t
##    tp = 1 - 0.5*t
##    a = (1 - t2).*tp
##    b = (t2 + t).*tp
##    c = (t2 - t).*tp/3
##    d = (t2 - 1).*t/6
##    J = matrix(1, 1,nc)
##    yi(range,:) = a(:,J) .* y(idx,:) + b(:,J) .* y(idx+1,:) ...
##        + c(:,J) .* y(idx-1,:) + d(:,J) .* y(idx+2,:)
##
##  elseif strcmp(method, 'spline') || strcmp(method, '*spline')
##    if (nx == 2) x = linspace(minx, maxx, ny); } #if
##    yi(range,:) = spline(x, y, xi)
##
##  else
##    error(["interp1 doesn't understand method '", method, "'"])
##  } #if
##  if transposed, yi=yi.'; } #if
##
##} #function

interp1 <- function(x, y, xi, 
    method = c('linear', 'nearest', 'pchip', 'cubic', 'spline'), 
    extrap = NA, ...) 
{
  method <- match.arg(method)
  lookup <- function(x, xi)
    approx(x, seq_along(x), xi, method = "constant", yleft = 0, yright = length(x))$y
  
  ny <- length(y)
  nx <- length(x)

  if (isTRUE(extrap == "extrap") || isTRUE(extrap)) {
    range <- seq_along(xi)
    yi <- numeric(length(xi))
  } else {
    range <- which(xi >= min(x) & xi <= max(x))
    yi <- rep.int(extrap, length(xi))
    xi <- xi[range]
  }

  yi[range] <- switch(method,
    "linear" = {
        if(nx < 2 || ny < 2)
            stop("interp1: table too short")
        if(nx == 2)
            idx <- rep(1, length(xi))
        else    
            idx <- lookup(x[2:(nx-1)], xi) + 1
                    # 2:(n-1) so that anything beyond the ends
                    # gets dumped into an interval
        ## use the endpoints of the interval to define a line
        dy <- diff(y)
        dx <- diff(x)
        s <- (xi - x[idx]) / dx[idx]
        s * dy[idx] + y[idx]
        #yi = approx(x, y, xi, ties = "ordered", ...)$y
    },
    "nearest" = {
        approx(c(x[1], x + c(diff(x)/2,0)), c(y[1], y), xi, method = "constant", f=1)$y
    }, 
    "spline" = {
        splinefun(x, y, ...)(xi)
    },
    "pchip" = {
        pchip(x, y, xi) 
    },
    "cubic" = {
        if (nx < 4 || ny < 4)
            stop("interp1: table too short")
        idx = lookup(x[3:(nx-2)], xi) + 1
    
        ## Construct cubic equations for each interval using divided
        ## differences (computation of c and d don't use divided differences
        ## but instead solve 2 equations for 2 unknowns). Perhaps
        ## reformulating this as a lagrange polynomial would be more efficient.
    
        i <- 1:(nx-3)
        dx <- diff(x)
        dx2 <- x[i+1]^2 - x[i]^2
        dx3 <- x[i+1]^3 - x[i]^3
        a <- diff(y, diff=3) / dx[i]^3 / 6
        b <- (diff(y[1:(nx-1)], diff=2) / dx[i]^2 - 6*a*x[i+1]) / 2
        c <- (diff(y[1:(nx-2)]) - a*dx3 - b*dx2) / dx[i]
        d <- y[i] - ((a*x[i] + b) * x[i] + c) * x[i]
        ((a[idx] * xi + b[idx]) * xi + c[idx]) * xi + d[idx]
    })
  return(yi)
}

#!demo
###xf=seq(0,10,length=500); yf = sin(2*pi*xf/5)
###xp=c(0,4,5,6,8,10); yp = sin(2*pi*xp/5)
###xp=c(0:4,6:10); yp = sin(2*pi*xp/5)
###xp=0:10; yp = sin(2*pi*xp/5)
###xp=c(0:1,3:10); yp = sin(2*pi*xp/5)
###lin=interp1(xp,yp,xf,'linear')
###spl=interp1(xp,yp,xf,'spline')
###cub=interp1(xp,yp,xf,'cubic')
###near=interp1(xp,yp,xf,'nearest')
###plot(xp,yp)
###lines(xf,lin,col="red")
###lines(xf,spl,col="green")
###lines(xf,cub,col="blue")
###lines(xf,near,col="purple")

#! plot(xf,yf,';original;',xf,near,';nearest;',xf,lin,';linear;',...
#!      xf,cub,';pchip;',xf,spl,';spline;',xp,yp,'*;;')
#! %--------------------------------------------------------
#! % confirm that interpolated function matches the original

#!demo
#! xf=0:0.05:10; yf = sin(2*pi*xf/5)
#! xp=0:10;      yp = sin(2*pi*xp/5)
#! lin=interp1(xp,yp,xf,'*linear')
#! spl=interp1(xp,yp,xf,'*spline')
#! cub=interp1(xp,yp,xf,'*cubic')
#! near=interp1(xp,yp,xf,'*nearest')
#! plot(xf,yf,';*original;',xf,near,';*nearest;',xf,lin,';*linear;',...
#!      xf,cub,';*cubic;',xf,spl,';*spline;',xp,yp,'*;;')
#! %--------------------------------------------------------
#! % confirm that interpolated function matches the original

#!shared xp, yp, xi, style
#! xp=0:5;      yp = sin(2*pi*xp/5)
#! xi = sort([-1, max(xp)*rand(1,6), max(xp)+1])

#!test style = 'nearest'
#!assert (interp1(xp, yp, [min(xp)-1, max(xp)+1]), [NaN, NaN])
#!assert (interp1(xp,yp,xp,style), yp, 100*eps)
#!assert (interp1(xp,yp,xp',style), yp', 100*eps)
#!assert (interp1(xp',yp',xp',style), yp', 100*eps)
#!assert (interp1(xp',yp',xp,style), yp, 100*eps)
#!assert (isempty(interp1(xp',yp',[],style)))
#!assert (isempty(interp1(xp,yp,[],style)))
#!assert (interp1(xp,[yp',yp'],xi(:),style),...
#!    [interp1(xp,yp,xi(:),style),interp1(xp,yp,xi(:),style)])
#!assert (interp1(xp,[yp',yp'],xi,style),
#!    interp1(xp,[yp',yp'],xi,['*',style]))

#!test style = 'linear'
#!assert (interp1(xp, yp, [-1, max(xp)+1]), [NaN, NaN])
#!assert (interp1(xp,yp,xp,style), yp, 100*eps)
#!assert (interp1(xp,yp,xp',style), yp', 100*eps)
#!assert (interp1(xp',yp',xp',style), yp', 100*eps)
#!assert (interp1(xp',yp',xp,style), yp, 100*eps)
#!assert (isempty(interp1(xp',yp',[],style)))
#!assert (isempty(interp1(xp,yp,[],style)))
#!assert (interp1(xp,[yp',yp'],xi(:),style),...
#!    [interp1(xp,yp,xi(:),style),interp1(xp,yp,xi(:),style)])
#!assert (interp1(xp,[yp',yp'],xi,style),
#!    interp1(xp,[yp',yp'],xi,['*',style]),100*eps)

#!test style = 'cubic'
#!assert (interp1(xp, yp, [-1, max(xp)+1]), [NaN, NaN])
#!assert (interp1(xp,yp,xp,style), yp, 100*eps)
#!assert (interp1(xp,yp,xp',style), yp', 100*eps)
#!assert (interp1(xp',yp',xp',style), yp', 100*eps)
#!assert (interp1(xp',yp',xp,style), yp, 100*eps)
#!assert (isempty(interp1(xp',yp',[],style)))
#!assert (isempty(interp1(xp,yp,[],style)))
#!assert (interp1(xp,[yp',yp'],xi(:),style),...
#!    [interp1(xp,yp,xi(:),style),interp1(xp,yp,xi(:),style)])
#!assert (interp1(xp,[yp',yp'],xi,style),
#!    interp1(xp,[yp',yp'],xi,['*',style]),1000*eps)

#!test style = 'spline'
#!assert (interp1(xp, yp, [-1, max(xp) + 1]), [NaN, NaN])
#!assert (interp1(xp,yp,xp,style), yp, 100*eps)
#!assert (interp1(xp,yp,xp',style), yp', 100*eps)
#!assert (interp1(xp',yp',xp',style), yp', 100*eps)
#!assert (interp1(xp',yp',xp,style), yp, 100*eps)
#!assert (isempty(interp1(xp',yp',[],style)))
#!assert (isempty(interp1(xp,yp,[],style)))
#!assert (interp1(xp,[yp',yp'],xi(:),style),...
#!    [interp1(xp,yp,xi(:),style),interp1(xp,yp,xi(:),style)])
#!assert (interp1(xp,[yp',yp'],xi,style),
#!    interp1(xp,[yp',yp'],xi,['*',style]),10*eps)

#!# test linear extrapolation
#!assert (interp1([1:5],[3:2:11],[0,6],'linear','extrap'), [1, 13], eps)
#!assert (interp1(xp, yp, [-1, max(xp)+1],'linear',5), [5, 5])

#!error interp1
#!error interp1(1:2,1:2,1,'bogus')

#!error interp1(1,1,1, 'nearest')
#!assert (interp1(1:2,1:2,1.4,'nearest'),1)
#!error interp1(1,1,1, 'linear')
#!assert (interp1(1:2,1:2,1.4,'linear'),1.4)
#!error interp1(1:3,1:3,1, 'cubic')
#!assert (interp1(1:4,1:4,1.4,'cubic'),1.4)
#!error interp1(1:2,1:2,1, 'spline')
#!assert (interp1(1:3,1:3,1.4,'spline'),1.4)

#!error interp1(1,1,1, '*nearest')
#!assert (interp1(1:2:4,1:2:4,1.4,'*nearest'),1)
#!error interp1(1,1,1, '*linear')
#!assert (interp1(1:2:4,1:2:4,[0,1,1.4,3,4],'*linear'),[NaN,1,1.4,3,NaN])
#!error interp1(1:3,1:3,1, '*cubic')
#!assert (interp1(1:2:8,1:2:8,1.4,'*cubic'),1.4)
#!error interp1(1:2,1:2,1, '*spline')
#!assert (interp1(1:2:6,1:2:6,1.4,'*spline'),1.4)
