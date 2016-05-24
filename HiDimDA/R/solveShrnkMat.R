### SolveShrnkMat.R  (2012-06-09)
###    
###
### Copyright 2012 A. Pedro Duarte Silva
###
### This file is part of the `HiDimDA' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

solve.ShrnkMat <- function(a,b=NULL,...)
{
  if (is.null(b))  {
	class(a) <- "ShrnkMatInv"
   	return(a) 
  }
  if (length(b)!=a$p) stop("Dimension of vector b is not compatible with dimensions of matrix a.\n")
  if (a$Trgt[1]=="Idntty") {
	ITrgtU <- a$U/a$Intst
	ITrgtb <- b/a$Intst
	if (a$q==1) InnerM <- 1./(a$D*(1-a$Intst)) + 1./a$Intst 
	else InnerM <- DMat(1./(a$D*(1-a$Intst))+rep(1./a$Intst,a$q))
  } 
  else {
	if (a$p==1) {
		ITrgt <- 1./(a$ITrgt*a$Intst)
		ITrgtb <- ITrgtb*b
		ITrgtU <- ITrgt*a$U
		tUITrgtb <- ITrgtU*b 
	}
	else {
		ITrgt <- solve(a$ITrgt)/a$Intst
		ITrgtb <- RightMult(ITrgt,b) 
		ITrgtU <- RightMult(ITrgt,a$U)
		tUITrgtb <- t(ITrgtU) %*% b 
	}
	if (a$q==1) {
		if (a$p==1) InnerM <- 1./(a$D*(1-a$Intst)) + ITrgt*a$U^2 
		else InnerM <- 1./(a$D*(1-a$Intst)) + t(a$U) %*% RightMult(ITrgt,a$U)
        }
	else InnerM <- diag(1./(a$D*(1-a$Intst))) + t(a$U) %*% RightMult(ITrgt,a$U)
  }
  if (a$q==1) return( drop(ITrgtb - ITrgtU*tUITrgtb/InnerM) )
  drop( ITrgtb - ITrgtU %*% solve(InnerM,tUITrgtb) ) 
}

solve.ShrnkMatInv <- function(a,b=NULL,...)
{
  if (is.null(b))  {
	class(a) <- "ShrnkMat"
   	return(a) 
  }
  if (length(b)!=a$p) stop("Dimension of vector b is not compatible with dimensions of matrix a.\n")
  if (a$q==1) UD <-  drop(a$U*a$D)
  else UD <- t(a$D*t(a$U))
  if (a$p==1) tUb <-  drop(a$U*b)
  else tUb <- t(a$U)%*%b
  if (a$Trgt[1]=="Idntty") {
	if (a$q==1) return( drop( a$Intst*b + (1-a$Intst)*UD*tUb ) )  
	return( drop( a$Intst*b + (1-a$Intst)*UD%*%tUb ) )  
  }
  if (a$q==1) {
	if (a$p==1) return( drop( a$Intst*a$Trgt*b + (1-a$Intst)*UD*tUb ) )  
	return( drop( a$Intst*RightMult(a$Trgt,b) + (1-a$Intst)*UD*tUb ) )
  }  
   drop( a$Intst*RightMult(a$Trgt,b) + (1-a$Intst)*(UD%*%tUb) )
} 


