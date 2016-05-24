### lw.dist.R                   
### Function to calculate conditional scores distributions using the Lord and
### Wingersky algorithm
###
### Requires a matrix It of item parameters estimates and the value Theta in
### to condition on
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

lw.dist<-function(It,Theta){
K<-dim(It)[2]
prov<-NULL
a<-It[2,1]
b<-It[1,1]
c<-It[3,1]
prov[1]<-1-Pr(Theta,b,a,c)
prov[2]<-Pr(Theta,b,a,c)
res<-prov
for (i in 2:K){
prov<-NULL
a<-It[2,i]
b<-It[1,i]
c<-It[3,i]
Pi<-Pr(Theta,b,a,c)
Qi<-1-Pi
prov[1]<-Qi*res[1]
for (j in 2:i) prov[j]<-Qi*res[j]+Pi*res[j-1]
prov[i+1]<-Pi*res[i]
res<-prov
}
return(res)}

