#  File degreenet/R/ghermite.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of California-Los Angeles
# Copyright 2007 The statnet Development Team
######################################################################
#
#  rmutil : A Library of Special Functions for Repeated Measurements
#  Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public Licence as published by
#  the Free Software Foundation; either version 2 of the Licence, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public Licence for more details.
#
#  You should have received a copy of the GNU General Public Licence
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     gauss.hermite(points, iterlim=10)
#
#  DESCRIPTION
#
#    Function to compute points and weights for Gauss-Hermite integration

gauss.hermite <- function(points, iterlim=10){
x <- w <- rep(0,points)
m <- (points+1)/2
for(i in 1:m){
	z <- if(i==1)sqrt(2*points+1)-2*(2*points+1)^(-1/6)
	else if(i==2)z-sqrt(points)/z
	else if(i==3||i==4)1.9*z-0.9*x[i-2]
	else 2.0*z-x[i-2]
	for(j in 1:iterlim){
		z1 <- z
		p <- hermite(points,z)
		z <- z1-p[1]/p[2]
		if(abs(z-z1)<=1e-15)break}
	if(j==iterlim)warning("iteration limit exceeded")
	x[points+1-i] <- -(x[i] <- z)
	w[i] <- w[points+1-i] <- 2/p[2]^2}
r <- cbind(x*sqrt(2),w/sum(w))
colnames(r) <- c("Points","Weights")
r}

# orthonormal Hermite polynomials
hermite <- function(points, z){
p1 <- 1/pi^0.4
p2 <- 0
for(j in 1:points){
	p3 <- p2
	p2 <- p1
	p1 <- z*sqrt(2.0/j)*p2-sqrt((j-1)/j)*p3}
pp <- sqrt(2*points)*p2
c(p1,pp)}
