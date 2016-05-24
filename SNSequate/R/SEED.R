### SEED.R
### Function to calculate the Standard Error of the Difference
### between two equating functions (Von Davier et al 2004 (p.80)
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

SEED<-function(eq1,eq2,...)
UseMethod("SEED")

SEED.default<-function(eq1,eq2,...)
{
	if(!is(eq1, "ker.eq") || !is(eq2, "ker.eq")){
		stop("eq1 and eq2 must be of class 'ker.eq'")}
else{
	if(eq1$design!="NEAT_CE" & eq2$design!="NEAT_CE"){
	RYx<-eq1$eqYx-eq2$eqYx
	sedVecYx<-eq1$sevecYx-eq2$sevecYx
	seedYx<-sqrt(apply(sedVecYx^2,1,sum))

	RXy<-eq1$eqXy-eq2$eqXy
	sedVecXy<-eq1$sevecXy-eq2$sevecXy
	seedXy<-sqrt(apply(sedVecXy^2,1,sum))

	res<- list(SEEDYx=seedYx,SEEDXy=seedXy,design=eq1$design)
	}
	else if(eq1$design=="NEAT_CE"){
	RYx<-eq1$eqCEYx-eq2$eqYx
	sedVecYx<-eq1$sevecYx-eq2$sevecYx
	seedYx<-sqrt(apply(sedVecYx^2,1,sum))
	res<- list(SEEDYx=seedYx,design=eq1$design)
	}
}
	class(res)<-"SEED"
res
}




print.SEED<-function(x,...)
{
	cat("\nStandard Error of Equating Difference between eq1 and eq2:\n")
	cat("\n")
	if(x$design!="NEAT_CE"){
	print(data.frame(SEEDYx=x$SEEDYx,SEEDXy=x$SEEDXy))
	}
	else{
	print(data.frame(SEEDYx=x$SEEDYx))
	}
}	



