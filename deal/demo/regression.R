## regression.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Mar 15 10:39:45 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Wed Jan 07 08:57:33 2004
## Update Count    : 16
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Susanne Gammelgaard Bottcher, Claus Dethlefsen
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################
n <- 1000

set.seed(109)

x <- seq(-2,2,length=n/2)
x2<- x^2

y1      <- rnorm(n/2,-5-x-x^2,.5)
y2      <- rnorm(n/2,+5+x+x^2,.5)
y       <- c(y1,y2)
A       <- factor(rep(c("A1","A2"),c(n/2,n/2)))
mypoly  <- data.frame(x,x2,y,A)
names(mypoly)[2] <- "x^2"

fit       <- network(mypoly)
fit.prior <- jointprior(fit)
res       <- nwfsort( getnetwork(networkfamily(mypoly,fit,fit.prior))  )

