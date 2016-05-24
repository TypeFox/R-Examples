### le.eq.R                   
### Function to implement the local equating method
###
### Requires the vector of observed scores "S.X" and the corresponding
### matrices of item parameter estimates for forms X and Y, as well as 
### the values on theta where to condition
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

le.eq<-function(S.X,It.X,It.Y,Theta) 
UseMethod("le.eq")

le.eq.default<-function(S.X,It.X,It.Y,Theta)
{
	###########################
	#Call parameters
	###########################
	cl<-match.call()
	

res<-tab<-NULL
U.T<-sort(unique(Theta))
for (i in 1:length(U.T)){
X.i<-S.X[Theta==U.T[i]]
U.X<-sort(unique(X.i))
pr<-NULL
for (j in 1:length(U.X)) {
pr<-c(U.T[i],U.X[j],LE(U.T[i],It.X,It.Y,U.X[j]))
tab<-rbind(tab,pr)
}
}
resu<-NULL
for (i in 1:length(S.X)) resu<-c(resu,tab[,3][tab[,1]==Theta[i] & tab[,2]==S.X[i]])
res<-list(call=cl,Theta=Theta,Obs.Sc=S.X,resu=as.numeric(resu))

class(res)<-"le.eq"
		   return(res)

}


print.le.eq<-function(x,...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\nEquated values under Local Equating:\n")	
	cat("\n")
	print(data.frame(Score=x$Obs.Sc,Th=x$Theta,eqYx=x$resu))
	cat("\n")
}	


