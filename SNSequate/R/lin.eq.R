### lin.eq.R                   
### Function that implments the linear method of equating
###
### Requires the vector of observed scores sx and sy and the value on the 
### scale to be equated
### 
### Copyright: Jorge Gonzalez, 2013.
### Last modification: 02-09-2013.
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

lin.eq<-function(sx,sy,scale) 
UseMethod("lin.eq")

lin.eq.default<-function(sx,sy,scale)
{
	###########################
	#Call parameters
	###########################
	cl<-match.call()

	mu.x<-mean(sx)
	mu.y<-mean(sy)
	sd.x<-sd(sx)
	sd.y<-sd(sy)
	resu<-(sd.y/sd.x)*(scale-mu.x)+mu.y
	res<-list(call=cl,scale=scale,resu=resu)
	class(res)<-"lin.eq"
	return(res)
}

print.lin.eq<-function(x,...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\nEquated values under linear Equating:\n")	
	cat("\n")
	print(data.frame(Score=x$scale,eqYx=x$resu))
	cat("\n")
}	