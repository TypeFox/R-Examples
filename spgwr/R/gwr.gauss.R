# Copyright 2001-2008 Roger Bivand

gwr.gauss <- function(dist2, bandwidth) {
	w <- exp((-dist2)/(bandwidth^2))
	w
}

gwr.Gauss <- function (dist2, bandwidth) {
     w <- exp((-0.5)*((dist2)/(bandwidth^2)))
     w
}

