### LE.R                   
### Performs equipercentile equating between the score distributions of 
### Forms X and Y, derived using the given value of ability
###
### Requires the value "Theta" of ability, the matrices "It.X" and "It.Y" of
### item parameters of both X and Y Forms, and the score "X" to equate.
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
### NOTE: The author thanks David Magis for providing a first version of 
###       this function


LE<-function(Theta,It.X,It.Y,X){
Dist.X<-lw.dist(It.X,Theta)
Dist.Y<-lw.dist(It.Y,Theta)
res<-EE(Dist.X,Dist.Y,X=X)
return(res)}
