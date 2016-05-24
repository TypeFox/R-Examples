###  Copyright (C) 2015 Berwin A. Turlach <Berwin.Turlach@gmail.com>
###
###  This program is free software; you can redistribute it and/or modify
###  it under the terms of the GNU General Public License as published by
###  the Free Software Foundation; either version 2 of the License, or
###  (at your option) any later version.
###
###  This program is distributed in the hope that it will be useful,
###  but WITHOUT ANY WARRANTY; without even the implied warranty of
###  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###  GNU General Public License for more details.
###
###  You should have received a copy of the GNU General Public License
###  along with this program; if not, write to the Free Software
###  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
###  USA.
###

ismonotone <- function(object, ...)
    UseMethod("ismonotone")


ismonotone.monpol <- function(object, a=-Inf, b=Inf, EPS=1e-6, ...){
    object <- coef(object)
    NextMethod("ismonotone")
}

ismonotone.default <- function(object, a=-Inf, b=Inf, EPS=1e-6, ...){
    ## need to get the derivitive
    deg <- length(object)-1
    deriv <- object[2:(deg+1)]*(1:deg)
    ## need to get the roots of the function
    roots <- polyroot(deriv)
    ## getting real roots only
    re.roots <- Re(roots)
    im.roots <- Im(roots)
    real.roots <- re.roots[abs(im.roots)<EPS]
    ## only real.roots between a and b
    real.roots <- real.roots[a < real.roots & real.roots < b] 
    ## checking multiplicity
    if(length(real.roots) == 0){
        TRUE
    }else{
        mltplcty <- rowSums(outer(real.roots, real.roots,
                                  function(x,y) abs(x-y)<EPS))
        all( mltplcty%%2 == 0)
    }
}
