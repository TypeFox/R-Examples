##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
##
##     Written by Karim Rahim and Wesley Burr.
##
##     This file is part of the multitaper package for R.
##     http://cran.r-project.org/web/packages/multitaper/index.html
## 
##     The multitaper package is free software: you can redistribute it and 
##     or modify it under the terms of the GNU General Public License as 
##     published by the Free Software Foundation, either version 2 of the 
##     License, or any later version.
##
##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
##
##     If you wish to report bugs please contact the author:
## 
##     Karim Rahim
##     karim.rahim@gmail.com


##############################################################
##
##  multitaperTrend 
##
##  Utility routine that computes multitaper-based linear
##  trend line. Has improved spectral properties over 
##  traditional least-squares. Returns intercept and slope.
##
##############################################################
multitaperTrend = function(xd, B, deltat, t.in) {

    N <- length(t.in)
    w <- B*deltat
    
    if(length(xd)!=N) { stop("Time array and data array not the same length!")} 
    if((B <= 0) || (B > 0.5)) { stop("B outside acceptable limits: 0 < B < 0.5.")}
    
    ttbar <- t.in - (t.in[N]+t.in[1])/2
    k <- floor(2*N*w -1)
    vt <- (dpss(N,k=k,nw=N*w))$v
    vk <- colSums(vt)
    
    ## solve for a
    subsel <- seq(1,k,by=2)
    vk <- colSums(vt)[subsel]
    xk <- colSums(xd*vt[,subsel])
    a <- sum(xk*vk) / sum(vk*vk)
    
    ## solve for b
    subsel <- seq(2,k,by=2)
    tvk <- colSums(ttbar*vt[,subsel])
    xk <- colSums(xd*vt[,subsel])
    b <- sum(tvk*xk)/sum(tvk*tvk)

    return(list(a,b,ttbar))
}
