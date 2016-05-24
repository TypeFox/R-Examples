### SigFqInv.R  (2012-05-21)
###    
###
### Copyright 2011 A. Pedro Duarte Silva
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

SigFqInv <- function(D,B,p,q,optres=NULL) solve(SigFq(D,B,p,q,optres))

is.SigFqInv <- function(x)  inherits(x,"SigFqInv")

as.matrix.SigFqInv <- function(x,...)  
{
   	if (is.null(x$B))  stop("Trying to convert a NULL SigFqInv object to a matrix.\n")
	if (x$p==1) return(matrix(x$D-x$B^2,1,1))
	diag(x$D)-x$B%*%t(x$B)  # return(diag(x$D)-x$B%*%t(x$B))
}

print.SigFqInv <- function(x,...)
{
	xi <- solve(x)
	cat(paste("Precision matrix for a ",x$q,"-factor model\n",sep=""))
	cat("\nLoadings:\n") ; print(xi$B)
	cat("\nSpecific Variances:\n",xi$D,"\n")
}

"+.SigFqInv" <- function(x,a)
{	
   if (is.SigFqInv(x)) {
   	if (is.null(x$B))  stop("Trying to add a NULL SigFq object\n")
   	if (is.SigFqInv(a)) {
		if (x$p!=a$p) stop("Argument dimensions do not match.\n")
   		if (x$p==1) return( matrix(x$D-x$B^2+a$D-a$B^2,1,1) )
   		return( diag(x$D)-x$B%*%t(x$B)+diag(a$D)-a$B%*%t(a$B) )
	}
   	if (!is.matrix(a) && length(a)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(a) && (x$p!=nrow(a) || x$p!=ncol(a)) ) stop("Argument dimensions do not match.\n")
        if (x$p==1)  return( matrix(x$D-x$B^2+a,1,1) )
        return( diag(x$D)-x$B%*%t(x$B)+a ) 
    }
    else {
   	if (is.null(a$B))  stop("Trying to add a NULL SigFq object\n")
   	if (!is.matrix(x) && length(x)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(x) && (a$p!=nrow(x) || a$p!=ncol(x)) ) stop("Argument dimensions do not match.\n")
        if (a$p==1)  return( matrix(a$D-a$B^2+x,1,1) )
        return( diag(a$D)-a$B%*%t(a$B)+x ) 
    }
}

"-.SigFqInv" <- function(x,a)
{	
   if (is.SigFqInv(x)) {
   	if (is.null(x$B))  stop("Trying to subtract from a NULL SigFq object\n")
   	if (is.SigFqInv(a)) {
		if (x$p!=a$p) stop("Argument dimensions do not match.\n")
   		if (x$p==1) return( matrix(x$D-x$B^2-a$D+a$B^2,1,1) )
   		return( diag(x$D)-x$B%*%t(x$B)-diag(a$D)+a$B%*%t(a$B) )
	}
   	if (!is.matrix(a) && length(a)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(a) && (x$p!=nrow(a) || x$p!=ncol(a)) ) stop("Argument dimensions do not match.\n")
        if (x$p==1)  return( matrix(x$D-x$B^2-a,1,1) )
        return( diag(x$D)-x$B%*%t(x$B)-a ) 
    }
    else {
   	if (is.null(a$B))  stop("Trying to subtract a NULL SigFq object\n")
   	if (!is.matrix(x) && length(x)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(x) && (a$p!=nrow(x) || a$p!=ncol(x)) ) stop("Argument dimensions do not match.\n")
        if (a$p==1)  return( matrix(x-a$D+a$B^2,1,1) )
        return( x-diag(a$D)+a$B%*%t(a$B) ) 
    }
}

"*.SigFqInv" <- function(x,a)
{
   if (is.SigFqInv(x)) {
   	if (is.null(x$B))  stop("Trying to multiply a NULL SigFq object\n")
   	if (is.SigFqInv(a)) {
		if (x$p!=a$p) stop("Argument dimensions do not match.\n")
   		if (x$p==1) return( matrix((x$D-x$B^2)*(a$D-a$B^2),1,1) )
   		return( (diag(x$D)-x$B%*%t(x$B))*(diag(a$D)-a$B%*%t(a$B)) )
	}
   	if (is.matrix(a)) {
		if (x$p!=nrow(a) || x$p!=ncol(a)) stop("Argument dimensions do not match.\n")
        	if (x$p==1)  return( matrix((x$D-x$B^2)*a,1,1) )
        	return( (diag(x$D)-x$B%*%t(x$B))*a )
	}
   	if (length(a)>1) stop("Arguments of wrong type.\n")
	result <- list(p=x$p,q=x$q,B=sqrt(a)*x$B,D=a*x$D,optres=NULL)
	class(result) <- "SigFqInv"
	return(result) 
    }
    else {
   	if (is.null(a$B))  stop("Trying to multiply a NULL SigFq object\n")
   	if (is.matrix(x)) {
		if (a$p!=nrow(x) || a$p!=ncol(x)) stop("Argument dimensions do not match.\n")
        	if (a$p==1)  return( matrix((a$D-a$B^2)*x,1,1) )
        	return( (diag(a$D)-a$B%*%t(a$B))*x )
	}
   	if (length(x)>1) stop("Arguments of wrong type.\n")
	result <- list(p=a$p,q=a$q,B=sqrt(x)*a$B,D=x*a$D,optres=NULL)
	class(result) <- "SigFqInv"
	return(result) 
    }
}

"/.SigFqInv" <- function(x,a)
{
   if (is.SigFqInv(x)) {
   	if (is.null(x$B))  stop("Trying to divide a NULL SigFq object\n")
   	if (is.SigFqInv(a)) {
		if (x$p!=a$p) stop("Argument dimensions do not match.\n")
   		if (x$p==1) return( matrix((x$D-x$B^2)/(a$D-a$B^2),1,1) )
   		return( (diag(x$D)-x$B%*%t(x$B))/(diag(a$D)-a$B%*%t(a$B)) )
	}
   	if (is.matrix(a)) {
		if (x$p!=nrow(a) || x$p!=ncol(a)) stop("Argument dimensions do not match.\n")
        	if (x$p==1)  return( matrix((x$D-x$B^2)/a,1,1) )
        	return( (diag(x$D)-x$B%*%t(x$B))/a )
	}
   	if (length(a)>1) stop("Arguments of wrong type.\n")
	result <- list(p=x$p,q=x$q,B=x$B/sqrt(a),D=x$D/a,optres=NULL)
	class(result) <- "SigFqInv"
	return(result) 
    }
    else {
   	if (is.null(a$B))  stop("Trying to divide from a NULL SigFq object\n")
   	if (!is.matrix(x) && length(x)>1) stop("Arguments of wrong type.\n")
   	if ( is.matrix(x) && (a$p!=nrow(x) || a$p!=ncol(x)) ) stop("Argument dimensions do not match.\n")
        if (a$p==1)  return( matrix(x/(a$D-a$B^2),1,1) )
        return( x/(diag(a$D)-a$B%*%t(a$B)) )
    }
}

LeftMult.SigFqInv <- function(x,a)
{
   if (is.SigFq(a) || is.SigFqInv(a)) {
	if (x$p!=a$p) stop("Argument dimensions do not match.\n")
	if (x$p>1) tmp1 <- diag(x$D*a$D)
	else tmp1 <- x$D*a$D
	tmp2 <- a$D * x$B 
	tmp3 <- x$D * a$B 
	tmp4 <- t(a$B) %*% x$B
	if (is.SigFq(a)) result <- tmp1 - tmp2%*%t(x$B) + a$B%*%t(tmp3) - a$B%*%tmp4%*%t(x$B) 
	if (is.SigFqInv(a)) result <- tmp1 - tmp2%*%t(x$B) - a$B%*%t(tmp3) + a$B%*%tmp4%*%t(x$B)
	return(result)
   }
   if (is.matrix(a))  { 
	if (x$p!=ncol(a)) stop("Argument dimensions do not match.\n")
   	t(t(a)*x$D) - (a%*%x$B) %*% t(x$B)  # return( t(t(a)*x$D) - (a%*%x$B) %*% t(x$B) )
   } 
   else  { 
	if (x$p!=length(a))  stop("Argument dimensions do not match.\n")
   	matrix(a*x$D,1,x$p) - (a%*%x$B) %*% t(x$B)  # return( matrix(a*x$D,1,x$p) - (a%*%x$B) %*% t(x$B) )
   }
}

RightMult.SigFqInv <- function(x,a)
{
   if (is.SigFq(a) || is.SigFqInv(a)) {
	if (x$p!=a$p) stop("Argument dimensions do not match.\n")
	if (x$p>1) tmp1 <- diag(x$D*a$D)
	else tmp1 <- x$D*a$D
	tmp2 <- x$D * a$B 
	tmp3 <- a$D * x$B 
	tmp4 <- t(x$B) %*% a$B
	if (is.SigFq(a)) result <- tmp1 + tmp2%*%t(a$B) - x$B%*%t(tmp3) - x$B%*%tmp4%*%t(a$B)  
	if (is.SigFqInv(a)) result <- tmp1 - tmp2%*%t(a$B) - x$B%*%t(tmp3) + x$B%*%tmp4%*%t(a$B)  
	return(result)
   }
   if (is.matrix(a))  { 
	if (x$p!=nrow(a)) stop("Argument dimensions do not match.\n")
   	a*x$D - x$B %*% (t(x$B)%*%a)  # return( a*x$D + (a%*%x$B) %*% t(x$B) )
   } 
   else  { 
	if (x$p!=length(a))  stop("Argument dimensions do not match.\n")
   	matrix(a*x$D,x$p,1) - x$B %*% (t(x$B)%*%a)  # return( matrix(a*x$D,x$p,1) - x$B %*% (t(x$B)%*%a) )
   }

}

