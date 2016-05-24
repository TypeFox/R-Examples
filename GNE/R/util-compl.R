#############################################################################
#   Copyright (c) 2012 Christophe Dutang                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################
### utility functions for complementarity functions in GNE
###
###         R functions
### 



#see page 857 or page xxx (30) of vol II of Facchinei & Pang (2003)


#Fischer-Burmeister
phiFB <- function(a, b) 
	sqrt(a^2+b^2) - (a+b)
GrAphiFB <- function(a, b) 
	ifelse(a == 0 & b == 0, -1/2, a / sqrt(a^2+b^2) - 1)
GrBphiFB <- function(a, b) 
	ifelse(a == 0 & b == 0, -1/2, b / sqrt(a^2+b^2) - 1)

#penalized Fischer-Burmeister
phipFB <- function(a, b, p) 
	sqrt(a^2+b^2) - (a+b) - p*pmax(a,0)*pmax(b,0)
GrAphipFB <- function(a, b, p) 
	ifelse(a == 0 & b == 0, -1/2, a / sqrt(a^2+b^2) - 1 - p*b*(a >= 0 & b >= 0))
GrBphipFB <- function(a, b, p) 
	ifelse(a == 0 & b == 0, -1/2, b / sqrt(a^2+b^2) - 1 - p*a*(a >= 0 & b >= 0))

#relaxed Fischer-Burmeister
phirFB <- function(a, b) 
	sqrt(a^2+b^2) - (a+b) - sqrt(sqrt(.Machine$double.eps))
GrAphirFB <- function(a, b) 
	ifelse(a == 0 & b == 0, -1/2, a / sqrt(a^2+b^2) - 1)
GrBphirFB <- function(a, b) 
	ifelse(a == 0 & b == 0, -1/2, b / sqrt(a^2+b^2) - 1)



#minimum
phiMin <- function(a, b) 
	min(a, b)
GrAphiMin <- function(a, b) 
	ifelse(a == 0 & b == 0, 1/2, 1*(a <= b))
GrBphiMin <- function(a, b) 
	ifelse(a == 0 & b == 0, 1/2, 1*(b <= a))


#Mangasarian
phiMan <- function(a, b, f, fprime)
	f(abs(a-b)) - f(a) - f(b)

GrAphiMan <- function(a, b, f, fprime)
	sign(a-b) * fprime(abs(a-b)) - fprime(a)

GrBphiMan <- function(a, b, f, fprime)
	sign(b-a) * fprime(abs(a-b)) - fprime(b)


#Luo-Tseng
phiLT <- function(a, b, q)
	(a^q + b^q)^(1/q) - (a+b)

GrAphiLT <- function(a, b, q)
	ifelse(a == 0 & b == 0,
		1 + (1/2)^((q-1)/q),
		1 - sign(a)*( abs(a) / (a^q + b^q)^(1/q) )^(q-1) )

GrBphiLT <- function(a, b, q)
	ifelse(a == 0 & b == 0,
		1 + (1/2)^((q-1)/q),
		1 - sign(b)*( abs(b) / (a^q + b^q)^(1/q) )^(q-1) )


#Kanzow-Kleinmichel
phiKK <- function(a, b, lambda)
	(sqrt( (a-b)^2 + 2*lambda*a*b ) - (a+b) ) / (2-lambda)

GrAphiKK <- function(a, b, lambda)
	ifelse( a == 0 & b == 0, 
		-1+sqrt(16*(2-lambda)*(lambda^2 - 2*lambda + 4))/8/(2-lambda),
		( (a+b*(lambda-1)) / sqrt( (a-b)^2 + 2*lambda*a*b ) - 1 ) / (2-lambda) )

GrBphiKK <- function(a, b, lambda)
	ifelse( a == 0 & b == 0, 
		   -1+sqrt(16*(2-lambda)*(lambda^2 - 2*lambda + 4))/8/(2-lambda),
		   ( (b+a*(lambda-1)) / sqrt( (a-b)^2 + 2*lambda*a*b ) - 1 ) / (2-lambda) )





#complementarity parameter class
compl.par <- function(type=c("FB", "pFB", "rFB", "Min", "Man", "LT", "KK"), 
	p, f, fprime, q, lambda)
{
	type <- match.arg(type, c("FB", "pFB", "rFB", "Min", "Man", "LT", "KK"))
	
	if(type == "FB")
		res <- list(type=type, fun=phiFB, grA=GrAphiFB, grB=GrBphiFB)
	else if(type == "pFB")
		res <- list(type=type, fun=phipFB, grA=GrAphipFB, grB=GrBphipFB, p=p)
	else if(type == "rFB")
		res <- list(type=type, fun=phirFB, grA=GrAphirFB, grB=GrBphirFB)
	else if(type == "Min")
		res <- list(type=type, fun=phiMin, grA=GrAphiMin, grB=GrBphiMin)
	else if(type == "Man")
		res <- list(type=type, fun=phiMan, grA=GrAphiMan, grB=GrBphiMan, f=f, fprime=fprime)
	else if(type == "LT")
		res <- list(type=type, fun=phiLT, grA=GrAphiLT, grB=GrBphiLT, q=q)
	else if(type == "KK")
		res <- list(type=type, fun=phiKK, grA=GrAphiKK, grB=GrBphiKK, lambda=lambda)
	else
		stop("internal error in compl.par.")
	
	class(res) <- "compl.par"
	res
}

#print function
print.compl.par <- function(x, ...)
{
	cat("Complementarity function:", x$type,"\n")
	print(x$fun)
	cat("with derivatives:\n")
	print(x$grA)
	print(x$grB)
	if(length(names(x)) > 4)
	{
		cat("with additional arguments:\n")
		print(x[!names(x) %in% c("fun", "grA", "grB", "type")])
	}
}


#summary function
summary.compl.par <- function(object, ...)
{
	structure(object, class = c("summary.compl.par", class(object)))	
}


#print function
print.summary.compl.par <- function(x, ...)
{
	cat("Complementarity function:", x$type,"\n")
	print(args(x$fun))
	cat("with derivatives:\n")
	print(args(x$grA))
	print(args(x$grB))
}
