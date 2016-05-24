### DMat.R  (2012-05-20)
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

DMat <- function(D) 
{
   if (is.matrix(D)) {
	if (nrow(D)!=ncol(D)) stop("D argument of DMat must be a vector or a square matrix\n.")
	if (nrow(D)>1)  D <- diag(D)
        else dim(D) <- NULL
   }
   if (!is.numeric(D)) stop("DMat argument of wrong type.\n")
   class(D) <- "DMat"
   D  # return(D) 
}

is.DMat <- function(x)  inherits(x,"DMat")

as.matrix.DMat <- function(x,...) 
{
	if (length(x)==1) return(matrix(x,1,1))
        diag(x)  # return(diag(x))
}

print.DMat <- function(x,...)
{
	cat("Diagonal Matrix:\n",x,"\n")
}

#drop.DMat <- function(x)
#{
#	if (length(x) > 1) stop("drop method is only available for one-dimensinoal DMat objects.\n")
#	as.numeric(x)
#}

"+.DMat" <- function(x,a)
{	
   if (is.DMat(x)) {
   	p <- length(x)
   	if (is.DMat(a)) {
		if (p!=length(a)) stop("Argument dimensions do not match.\n")
		return(DMat(as.numeric(x)+as.numeric(a)))
   	}
   	if (!is.matrix(a) && length(a)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(a) && (p!=nrow(a) || p!=ncol(a)) ) stop("Argument dimensions do not match.\n")
   	if (p==1) return(matrix(as.numeric(x)+as.numeric(a),1,1))
   	return( as.matrix(x)+a )
   }
   else {
   	p <- length(a)
   	if (p==1) return(matrix(as.numeric(a)+as.numeric(x),1,1))
   	if (!is.matrix(x) && length(x)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(x) && (p!=nrow(x) || p!=ncol(x)) ) stop("Argument dimensions do not match.\n")
   	return( as.matrix(a)+x )
   }
}

"-.DMat" <- function(x,a)
{	
   if (is.DMat(x)) {
   	p <- length(x)
   	if (is.DMat(a)) {
		if (p!=length(a)) stop("Argument dimensions do not match.\n")
		return(DMat(as.numeric(x)-as.numeric(a)))
   	}
   	if (p==1) return(matrix(as.numeric(x)-as.numeric(a),1,1))
   	if (!is.matrix(a) && length(a)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(a) && (p!=nrow(a) || p!=ncol(a)) ) stop("Argument dimensions do not match.\n")
   	return( as.matrix(x)-a )
   }
   else {
   	p <- length(a)
   	if (p==1) return(matrix(as.numeric(x)-as.numeric(a),1,1))
   	if (!is.matrix(x) && length(x)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(x) && (p!=nrow(x) || p!=ncol(x)) ) stop("Argument dimensions do not match.\n")
   	return( x-as.matrix(a) )
   }
}

"*.DMat" <- function(x,a)
{
   if (is.DMat(x)) {
   	p <- length(x)
   	if (is.DMat(a)) {
		if (p!=length(a)) stop("Argument dimensions do not match.\n")
		return(DMat(as.numeric(x)*as.numeric(a)))
   	}
   	if (p==1) return(DMat(as.numeric(x)*as.numeric(a)))
   	if (!is.matrix(a) && length(a)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(a) && (p!=nrow(a) || p!=ncol(a)) ) stop("Argument dimensions do not match.\n")
   	if (!is.matrix(a)) return( DMat(as.numeric(x)*a) )
   	return( DMat(as.numeric(x)*diag(a)) )
   }
   else {
   	p <- length(a)
   	if (p==1) return(DMat(as.numeric(a)*as.numeric(x)))
   	if (!is.matrix(x) && length(x)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(x) && (p!=nrow(x) || p!=ncol(x)) ) stop("Argument dimensions do not match.\n")
   	if (!is.matrix(x)) return( DMat(as.numeric(a)*x) )
	return( DMat(as.numeric(a)*diag(x)) )
   }
}

"/.DMat" <- function(x,a)
{
   if (is.DMat(x)) {
   	p <- length(x)
   	if (is.DMat(a)) {
		if (p!=length(a)) stop("Argument dimensions do not match.\n")
		return(DMat(as.numeric(x)/as.numeric(a)))
   	}
   	if (!is.matrix(a) && length(a)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(a) && (p!=nrow(a) || p!=ncol(a)) ) stop("Argument dimensions do not match.\n")
   	if (p==1) return(DMat(as.numeric(x)/as.numeric(a)))
   	if (!is.matrix(a)) return( DMat(as.numeric(x)/a) )
   	return( DMat(as.numeric(x)/diag(a)) )
   }
   else {
   	p <- length(a)
   	if (p==1) return(DMat(as.numeric(x)/as.numeric(a)))
   	if (!is.matrix(x) && length(x)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(x) && (p!=nrow(x) || p!=ncol(x)) ) stop("Argument dimensions do not match.\n")
   	if (!is.matrix(x)) return( DMat(x/as.numeric(a)) )
	return( DMat(diag(x)/as.numeric(a)) )
   }
}

LeftMult.DMat <- function(x,a)
{
   p <- length(x)
   if (is.DMat(a)) {
	if (p!=length(a)) stop("Argument dimensions do not match.\n")
	return(DMat(as.numeric(x)*as.numeric(a)))
   }
   if (is.matrix(a))  { 
	if (p!=ncol(a)) stop("Argument dimensions do not match.\n")
   	sweep(a,2,as.numeric(x),FUN="*")  # return( sweep(a,2,as.numeric(x),FUN="*")) )
   } 
   else  { 
	if (p!=length(a))  stop("Argument dimensions do not match.\n")
   	matrix(as.numeric(x)*a,1,p)  # return( matrix(as.numeric(x)*a,1,p) )
   }
}

RightMult.DMat <- function(x,a)
{
   p <- length(x)
   if (is.DMat(a)) {
	if (p!=length(a)) stop("Argument dimensions do not match.\n")
	return(DMat(as.numeric(x)*as.numeric(a)))
   }
   if (is.matrix(a))  { 
	if (p!=nrow(a)) stop("Argument dimensions do not match.\n")
   	sweep(a,1,as.numeric(x),FUN="*")  # return( sweep(a,1,as.numeric(x),FUN="*")) )
   } 
   else  { 
	if (p!=length(a))  stop("Argument dimensions do not match.\n")
   	matrix(as.numeric(x)*a,p,1)  # return( matrix(as.numeric(x)*a,p,1) )
   }
}

solve.DMat <- function(a,b=NULL,...)
{
  if (is.null(b)) {
	result <- 1./as.numeric(a)
	class(result) <- "DMat"
  }
  else result <- b/as.numeric(a)
  result # return(result)
}

