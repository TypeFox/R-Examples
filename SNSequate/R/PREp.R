### PREp.R
### Function to calculate the Percent Relative Error
### according to the description in Von Davier et al 2004 (p.66)
###
### Copyright: Jorge Gonzalez, 2012.
### Last modification: 25-05-2012.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### Author contact information:
###
###      Jorge Gonzalez B.
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3545467  URL  : http://www.mat.puc.cl/~jgonzale
###      Fax  : +56-2-3547729  Email: jgonzale@mat.puc.cl
###
PREp<-function(eq,p)
UseMethod("PREp")

PREp.default<-function(eq,p){
  if (!is(eq, "ker.eq")) {
    stop("eq must be of class 'ker.eq'")
  }

  if (eq$design != "NEAT_CE") {
    preYx <- c()
    preXy <- c()
    
    for (i in 1:p) {
      preYx[i] <- 100 * (sum(eq$eqYx^i * eq$rj) - sum(eq$score^i * eq$sk)) / sum(eq$score^i * eq$sk)
      preXy[i] <- 100 * (sum(eq$eqXy^i * eq$sk) - sum(eq$score^i * eq$rj)) / sum(eq$score^i * eq$rj)
    }
    res <- list(Moments = 1:p,preYx = preYx,preXy = preXy,design = eq$design)
  }
  else {
    preYx <- c()
    for (i in 1:p) {
      preYx[i] <- 100 * (sum(eq$eqYx ^ i * eq$rj) - sum(eq$score^i * eq$sk)) / sum(eq$score^i * eq$sk)
    }
    res <- list(Moments = 1:p,preYx = preYx,design = eq$design)
  }
  
  class(res) <- "PREp"
  res
}



print.PREp<-function(x,...)
{
	cat("\nPercent Relative Error:\n")
	cat("\n")
	if(x$design!="NEAT_CE"){
	print(data.frame(Moments=x$Moments,'X_to_Y'=x$preYx,'Y_to_X'=x$preXy))
	}
}	




