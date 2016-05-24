### cdf.Dist.R                   
### Function to create CDF of scores to be used in equipercentile equating
###
### Requires the vector of observed scores "S.X" and the corresponding
### frequencies of each observed score "Dist.X"
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


cdf.Dist<-function(Dist.X,S.X){
cu<-cumsum(Dist.X)
prov.seq<-c(S.X[1]-0.5,S.X+0.5)
prov.dist<-c(0,cu)
res.seq<-seq(min(prov.seq),max(prov.seq),1)
res.dist<-NULL
for (i in 1:length(res.seq)){
test<-F
for (j in 1:length(prov.seq)){
if (prov.seq[j]==res.seq[i]) test<-T}
if (test==T) res.dist<-c(res.dist,prov.dist[prov.seq==res.seq[i]])
else res.dist<-c(res.dist,max(prov.dist[prov.seq<=res.seq[i]]))
}
res<-list(res.seq=res.seq,res.dist=res.dist)
return(res)
}
