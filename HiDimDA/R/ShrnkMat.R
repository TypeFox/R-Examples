### ShrnkMat.R  (2012-06-04)
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

ShrnkMat <- function(U,D,p,q,Intst,Trgt="Idntty") 
{
   if (!is.matrix(U)) 
	if (q==1 && length(U)==p) U <- matrix(U,p,1) 
	else stop("U argument of ShrnkMat must be a matrix.\n")   
   if (length(D)!=q || nrow(U)!=p || ncol(U)!=q) stop("Wrong argument dimensions.\n")
   if (is.matrix(Trgt) && !isTRUE(all.equal(Trgt,t(Trgt)))) 
	stop("ShrnkMat only accpets symmetric matrix targets\n")
   if (p<q) stop("Rank of ShrnkMat cannot be higher than its dimension.\n")
   result <- list(p=p,q=q,Trgt=Trgt,U=U,D=D,Intst=Intst)
   class(result) <- "ShrnkMat"
   result  # return(result) 
}

is.ShrnkMat <- function(x)  inherits(x,"ShrnkMat")

as.matrix.ShrnkMat<- function(x,...) 
{
   if (x$Trgt[1]=="Idntty")  {
	if (x$p==1) return(matrix((1.-x$Intst)*x$D*x$U^2+x$Intst,1,1))
	if (x$q==1) return( (1.-x$Intst)*(x$D*x$U%*%t(x$U)) + x$Intst*diag(x$p) )
	return( (1.-x$Intst) * (x$U%*%diag(x$D)%*%t(x$U)) + x$Intst*diag(x$p) )
   }  
   if (x$p==1) return(matrix((1.-x$Intst)*x$D*x$U^2+x$Intst*x$Trgt,1,1))
   if (x$q==1) return( (1.-x$Intst)*(x$D*x$U%*%t(x$U)) + x$Intst*as.matrix(x$Trgt) )
   (1.-x$Intst) * (x$U%*%diag(x$D)%*%t(x$U)) + x$Intst*as.matrix(x$Trgt)  
}

print.ShrnkMat <- function(x,...)
{
	cat("Shrunken symmetric Matrix estimator:\n")
	cat("Matrix dimension:",x$p,"\n")
	cat("Original Matrix rank:",x$q,"\n")
	cat("Original Matrix eigenvalues:\n") ; print(x$D) 
	cat("Original Matrix eigenvectors:\n") ; print(x$U)
	if (x$Trgt[1]=="Idntty") cat(paste("Target: I_",x$p,"\n",sep="")) 
	else { cat("Target:\n") ; print(x$Trgt) } 
	cat("Target Intensity:",x$Intst,"\n") 
}

"+.ShrnkMat" <- function(x,a)
{	
   if (is.ShrnkMat(x)) {
   	if (is.ShrnkMat(a)) {
		if (x$p!=a$p) stop("Argument dimensions do not match.\n")
		return(as.matrix(x)+as.matrix(a))
   	}
   	if (!is.matrix(a) && length(a)>1) stop("Arguments of wrong type.\n")
        return( as.matrix(x)+a )
    }
    else {
   	if (!is.matrix(x) && length(x)>1) stop("Arguments of wrong type.\n")
        return( x+as.matrix(a) )
    }
}

"-.ShrnkMat" <- function(x,a)
{	
   if (is.ShrnkMat(x)) {
   	if (is.ShrnkMat(a)) {
		if (x$p!=a$p) stop("Argument dimensions do not match.\n")
		return(as.matrix(x)-as.matrix(a))
   	}
   	if (!is.matrix(a) && length(a)>1) stop("Arguments of wrong type.\n")
        return( as.matrix(x)-a )
    }
    else {
   	if (!is.matrix(x) && length(x)>1) stop("Arguments of wrong type.\n")
        return( x-as.matrix(a) )
    }
}

"*.ShrnkMat" <- function(x,a)
{
   if (is.ShrnkMat(x)) {
   	if (is.ShrnkMat(a)) {
		if (x$p!=a$p) stop("Argument dimensions do not match.\n")
   		return(as.matrix(x)*as.matrix(a))
	}
   	if (is.matrix(a)) {
		if (x$p!=nrow(a) || x$p!=ncol(a)) stop("Argument dimensions do not match.\n")
   		return(as.matrix(x)*a)
	}
   	if (length(a)>1) stop("Arguments of wrong type.\n")
	if (x$Trgt[1]=="Idntty") result <- list(p=x$p,q=x$q,Trgt=a*diag(x$p),U=x$U,D=a*x$D,Intst=x$Intst)
	else result <- list(p=x$p,q=x$q,Trgt=a*x$Trgt,U=x$U,D=a*x$D,Intst=x$Intst)
	class(result) <- "ShrnkMat"
	return(result) 
    }
    else {
   	if (is.matrix(x)) {
		if (a$p!=nrow(x) || a$p!=ncol(x)) stop("Argument dimensions do not match.\n")
   		return(x*as.matrix(a))
	}
   	if (length(x)>1) stop("Arguments of wrong type.\n")
	if (a$Trgt[1]=="Idntty") result <- list(p=a$p,q=a$q,Trgt=x*diag(a$p),U=a$U,D=x*a$D,Intst=a$Intst)
	else result <- list(p=a$p,q=a$q,Trgt=x*a$Trgt,U=a$U,D=x*a$D,Intst=a$Intst)
	class(result) <- "ShrnkMat"
	return(result) 
    }
}

"/.ShrnkMat" <- function(x,a)
{
   if (is.ShrnkMat(x)) {
   	if (is.ShrnkMat(a)) {
		if (x$p!=a$p) stop("Argument dimensions do not match.\n")
   		return(as.matrix(x)/as.matrix(a))
	}
   	if (is.matrix(a)) {
		if (x$p!=nrow(a) || x$p!=ncol(a)) stop("Argument dimensions do not match.\n")
   		return(as.matrix(x)/a)
	}
   	if (length(a)>1) stop("Arguments of wrong type.\n")
	if (x$Trgt[1]=="Idntty") result <- list(p=x$p,q=x$q,Trgt=diag(x$p)/a,U=x$U,D=x$D/a,Intst=x$Intst)
	else result <- list(p=x$p,q=x$q,Trgt=x$Trgt/a,U=x$U,D=x$D/a,Intst=x$Intst)
	class(result) <- "ShrnkMat"
	return(result) 
    }
    else {
   	if (is.matrix(x)) {
		if (a$p!=nrow(x) || a$p!=ncol(x)) stop("Argument dimensions do not match.\n")
   		return(x/as.matrix(a))
	}
   	if (length(x)>1) stop("Arguments of wrong type.\n")
	return(x/as.matrix(a))
    }
}

LeftMult.ShrnkMat<- function(x,a)
{
   if (is.ShrnkMat(a) || is.ShrnkMatInv(a)) {
	if (x$p!=a$p) stop("Argument dimensions do not match.\n")
	A <- as.matrix(a)
	if (x$Trgt[1]=="Idntty") return( (1.-x$Intst)*(A%*%x$U)%*%sweep(t(x$U),1,x$D,FUN="*") + x$Intst*A )
	else return( (1.-x$Intst)*(A%*%x$U)%*%sweep(t(x$U),1,x$D,FUN="*") + x$Intst*LeftMult(x$Trgt,A) )
   }
   if (is.matrix(a))  { 
	if (x$p!=ncol(a)) stop("Argument dimensions do not match.\n")
	if (x$Trgt[1]=="Idntty") return( (1.-x$Intst)*(a%*%x$U)%*%sweep(t(x$U),1,x$D,FUN="*") + x$Intst*a )
	else return( (1.-x$Intst)*(a%*%x$U)%*%sweep(t(x$U),1,x$D,FUN="*") + x$Intst*LeftMult(x$Trgt,a) )
   } 
   if (x$p!=length(a))  stop("Argument dimensions do not match.\n")
   if (x$p!=1) dim(a) <- c(1,x$p)
   if (x$Trgt[1]=="Idntty")  {
	if (x$p==1) return( matrix((1.-x$Intst)*a*x$D*x$U^2+x$Intst*a,length(a),1) )
	if (x$q==1)  return( (1.-x$Intst) * (x$D*(a%*%x$U)%*%t(x$U)) + x$Intst*a )
	return( (1.-x$Intst) * ((a%*%x$U)%*%diag(x$D)%*%t(x$U)) + x$Intst*a )
   } 
   if (x$p==1) return(matrix((1.-x$Intst)*a*drop(x$D)*drop(x$U)^2+x$Intst*a*drop(x$Trgt),length(a),1) )
   if (x$q==1)  return( (1.-x$Intst) * (x$D*(a%*%x$U)%*%t(x$U)) + x$Intst*a%*%as.matrix(x$Trgt) )  
   (1.-x$Intst) * ((a%*%x$U)%*%diag(x$D)%*%t(x$U)) + x$Intst*a%*%as.matrix(x$Trgt)  
}

RightMult.ShrnkMat<- function(x,a)
{
   if (is.ShrnkMat(a) || is.ShrnkMatInv(a)) {
	if (x$p!=a$p) stop("Argument dimensions do not match.\n")
	A <- as.matrix(a)
	if (x$Trgt[1]=="Idntty") return( (1.-x$Intst)*sweep(x$U,2,x$D,FUN="*")%*%(t(x$U)%*%A) + x$Intst*A )
	else return( (1.-x$Intst)*sweep(x$U,2,x$D,FUN="*")%*%(t(x$U)%*%A) + x$Intst*RightMult(x$Trgt,A) )
   }
   if (is.matrix(a))  { 
	if (x$p!=nrow(a)) stop("Argument dimensions do not match.\n")
	if (x$Trgt[1]=="Idntty") return( (1.-x$Intst)*sweep(x$U,2,x$D,FUN="*")%*%(t(x$U)%*%a) + x$Intst*a )
	else return( (1.-x$Intst)*sweep(x$U,2,x$D,FUN="*")%*%(t(x$U)%*%a) + x$Intst*RightMult(x$Trgt,a) )
   } 
   if (x$p!=length(a))  stop("Argument dimensions do not match.\n")
   if (x$p!=1) dim(a) <- c(x$p,1)
   if (x$Trgt[1]=="Idntty")  {
	if (x$p==1) return(matrix((1.-x$Intst)*a*drop(x$D)*drop(x$U)^2+x$Intst*a,1,length(a)) )
	if (x$q==1)  return( (1.-x$Intst) * (x$D*x$U%*%t(x$U)%*%a) + x$Intst*a )
	return( (1.-x$Intst) * (x$U%*%diag(x$D)%*%t(x$U)%*%a) + x$Intst*a )
   } 
   if (x$p==1) return(matrix((1.-x$Intst)*a*x$D*x$U^2+x$Intst*a*x$Trgt,1,length(a)) )
   if (x$q==1)  return( (1.-x$Intst) * (x$D*x$U%*%t(x$U)%*%a) + x$Intst*as.matrix(x$Trgt)%*%a )  
   (1.-x$Intst) * (x$U%*%diag(x$D)%*%t(x$U)%*%a) + x$Intst*as.matrix(x$Trgt)%*%a  
}


