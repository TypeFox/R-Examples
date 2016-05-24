### ShrnkMatInv.R  (2012-06-08)
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

ShrnkMatInv<- function(U,D,p,q,Intst,Trgt="Idntty") 
{
   if (!is.matrix(U)) 
	if (q==1 && length(U)==p) U <- matrix(U,p,1) 
	else stop("U argument of ShrnkMat must be a matrix.\n")   
   if (length(D)!=p || nrow(U)!=p || ncol(U)!=q) stop("Wrong argument dimensions.\n")
   if (is.matrix(Trgt) && !isTRUE(all.equal(Trgt,t(Trgt)))) 
	stop("ShrnkMatInv only accepts symmetric matrix targets.\n")
   if (p<q) stop("Rank of ShrnkMatInv cannot be higher than its dimension.\n")
   result <- list(p=p,q=q,Trgt=Trgt,U=U,D=D,Intst=Intst)
   class(result) <- "ShrnkMatInv"
   result  # return(result) 
}

is.ShrnkMatInv<- function(x)  inherits(x,"ShrnkMatInv")

as.matrix.ShrnkMatInv<- function(x,...) 
{
  if (x$Trgt[1]=="Idntty") {
	ITrgt <- diag(x$p)/x$Intst
	ITrgtU <- x$U/x$Intst
	if (x$p==1) InnerM <- 1./(x$D*(1-x$Intst)) + (x$U^2)/x$Intst 
	else if (x$q==1) InnerM <- 1./(x$D*(1-x$Intst)) + 1./x$Intst 
	else InnerM <- DMat(1./(x$D*(1-x$Intst))+rep(1./x$Intst,x$q))
  } 
  else {
	if (x$p==1) {
		ITrgt <- 1./(x$Trgt*x$Intst)
		ITrgtU <- ITrgt*x$U
	}
	else {
		ITrgt <- solve(x$Trgt)/x$Intst
		ITrgtU <- RightMult(ITrgt,x$U)
	}
	if (x$q==1) {
		if (x$p==1) InnerM <- 1./(x$D*(1-x$Intst)) + ITrgt*x$U^2 
		else InnerM <- 1./(x$D*(1-x$Intst)) + drop(t(x$U)%*%ITrgtU)
        }
	else InnerM <- diag(1./(x$D*(1-x$Intst))) + t(x$U) %*% ITrgtU
  }
  if (x$q==1) {
	if (x$p==1) return( ITrgt - ITrgtU^2/InnerM )
	return( ITrgt - ITrgtU%*%t(ITrgtU)/InnerM )
  }
  ITrgt - ITrgtU %*% solve(InnerM,t(ITrgtU))
}

print.ShrnkMatInv<- function(x,...)
{
	cat("Inverse Shrunken symmetric Matrix estimator:\n")
	cat("Matrix dimension:",x$p,"\n")
	cat("Original Matrix rank:",x$q,"\n")
	cat("Original Matrix eigenvalues:\n") ; print(x$D) 
	cat("Original Matrix eigenvectors:\n") ; print(x$U)
	if (x$Trgt[1]=="Idntty") cat(paste("Target: I_",x$p,"\n",sep="")) 
	else { cat("Target:\n") ; print(x$Trgt) } 
	cat("Target Intensity:",x$Intst,"\n") 
}

"+.ShrnkMatInv" <- function(x,a)
{	
   if (is.ShrnkMatInv(x)) {
   	if (is.ShrnkMatInv(a)) {
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

"-.ShrnkMatInv" <- function(x,a)
{	
   if (is.ShrnkMatInv(x)) {
   	if (is.ShrnkMatInv(a)) {
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

"*.ShrnkMatInv" <- function(x,a)
{
   if (is.ShrnkMatInv(x)) {
   	if (is.ShrnkMatInv(a)) {
		if (x$p!=a$p) stop("Argument dimensions do not match.\n")
   		return(as.matrix(x)*as.matrix(a))
	}
   	if (is.matrix(a)) {
		if (x$p!=nrow(a) || x$p!=ncol(a)) stop("Argument dimensions do not match.\n")
   		return(as.matrix(x)*a)
	}
   	if (length(a)>1) stop("Arguments of wrong type.\n")
	if (x$Trgt[1]=="Idntty") x$Trgt <- diag(x$p)
	result <- list(p=x$p,q=x$q,Trgt=x$Trgt/a,U=x$U,D=x$D/a,Intst=x$Intst)
	class(result) <- "ShrnkMatInv"
	return(result) 
    }
    else {
   	if (is.matrix(x)) {
		if (a$p!=nrow(x) || a$p!=ncol(x)) stop("Argument dimensions do not match.\n")
   		return(x*as.matrix(a))
	}
   	if (length(x)>1) stop("Arguments of wrong type.\n")
	if (a$Trgt[1]=="Idntty") return(x*as.matrix(a))
	result <- list(p=a$p,q=a$q,Trgt=a$Trgt/x,U=a$U,D=a$D/x,Intst=a$Intst)
	class(result) <- "ShrnkMatInv"
	return(result) 
    }
}

"/.ShrnkMatInv" <- function(x,a)
{
   if (is.ShrnkMatInv(x)) {
   	if (is.ShrnkMatInv(a)) {
		if (x$p!=a$p) stop("Argument dimensions do not match.\n")
   		return(as.matrix(x)/as.matrix(a))
	}
   	if (is.matrix(a)) {
		if (x$p!=nrow(a) || x$p!=ncol(a)) stop("Argument dimensions do not match.\n")
   		return(as.matrix(x)/a)
	}
   	if (length(a)>1) stop("Arguments of wrong type.\n")
	if (x$Trgt[1]=="Idntty") x$Trgt <- diag(x$p)
	result <- list(p=x$p,q=x$q,Trgt=x$Trgt*a,U=x$U,D=x$D*a,Intst=x$Intst)
	class(result) <- "ShrnkMatInv"
	return(result) 
    }
    else {
   	if (is.matrix(x) && (a$p!=nrow(x) || a$p!=ncol(x)) ) stop("Argument dimensions do not match.\n")
   	if (!is.matrix(x) && length(x)>1) stop("Arguments of wrong type.\n")
  	return(x/as.matrix(a))
    }
}

LeftMult.ShrnkMatInv<- function(x,a)
{
  if (x$Trgt[1]=="Idntty") {
	ITrgt <- diag(x$p)/x$Intst
	ITrgtU <- x$U/x$Intst
	if (x$q==1) InnerM <- 1./(drop(x$D)*(1-x$Intst)) + 1./x$Intst 
	else InnerM <- DMat(1./(x$D*(1-x$Intst))+rep(1./x$Intst,x$q))
  } 
  else {
	if (x$p==1) {
		ITrgt <- 1./(drop(as.matrix(x$Trgt))*x$Intst)
		ITrgtU <- ITrgt*drop(x$U)
	}
	else {
		ITrgt <- solve(x$Trgt)/x$Intst
		ITrgtU <- RightMult(ITrgt,x$U)
	}
	if (x$q==1) {
		if (x$p==1) InnerM <- 1./(drop(x$D)*(1-x$Intst)) + ITrgt*drop(x$U)^2 
		else InnerM <- 1./(x$D*(1-x$Intst)) + drop(t(x$U)%*%RightMult(ITrgt,x$U))
        }
	else InnerM <- diag(1./(x$D*(1-x$Intst))) + t(x$U) %*% RightMult(ITrgt,x$U)
  }
  if (is.ShrnkMatInv(a) || is.ShrnkMat(a)) {
	if (x$p!=a$p) stop("Argument dimensions do not match.\n")
  	if (x$p==1) {
		aITrgt <- ITrgt*as.matrix(a)
		aITrgtU <- ITrgtU*as.matrix(a)
	} 
 	else aITrgt <- {
		LeftMult(ITrgt,as.matrix(a)) 
		aITrgtU <- as.matrix(a) %*% ITrgtU
	} 
  }
  else {
	if (is.matrix(a)) {
		if (x$p!=ncol(a)) stop("Argument dimensions do not match.\n")
        }
        else {
		if (x$p!=length(a))  stop("Argument dimensions do not match.\n")
		dim(a) <- c(1,x$p)
	}
  	if (x$p==1) {
		aITrgt <- matrix(ITrgt*drop(a),length(a),1) 
		aITrgtU <- matrix(ITrgtU*drop(a),length(a),1)  
	}
  	else {
		aITrgt <- LeftMult(ITrgt,a) 
		aITrgtU <- a %*% ITrgtU
	} 
   }
   if (x$q==1) {
	if (x$p==1) return( aITrgt - aITrgtU*drop(ITrgtU)/InnerM )
	return( aITrgt - aITrgtU%*%t(ITrgtU)/InnerM )
   }
   if (x$p==1) return(aITrgt - aITrgtU %*% ITrgtU/InnerM)
   aITrgt - aITrgtU %*% solve(InnerM,t(ITrgtU))
}

RightMult.ShrnkMatInv<- function(x,a)
{
  if (x$Trgt[1]=="Idntty") {
	ITrgt <- diag(x$p)/x$Intst
	ITrgtU <- x$U/x$Intst
	if (x$q==1) InnerM <- 1./(drop(x$D)*(1-x$Intst)) + 1./x$Intst
	else InnerM <- DMat(1./(x$D*(1-x$Intst))+rep(1./x$Intst,x$q))
  } 
  else {
	if (x$p==1) {
		ITrgt <- 1./(drop(as.matrix(x$Trgt))*x$Intst)
		ITrgtU <- ITrgt*drop(x$U)
	}
	else {
		ITrgt <- solve(x$Trgt)/x$Intst
		ITrgtU <- RightMult(ITrgt,x$U)
	}
	if (x$q==1) {
		if (x$p==1) InnerM <- 1./(drop(x$D)*(1-x$Intst)) + ITrgt*drop(x$U)^2 
		else InnerM <- 1./(x$D*(1-x$Intst)) + drop(t(x$U)%*%RightMult(ITrgt,x$U))
        }
	else InnerM <- diag(1./(x$D*(1-x$Intst))) + t(x$U) %*% RightMult(ITrgt,x$U)
  }
  if (is.ShrnkMatInv(a) || is.ShrnkMat(a)) {
	if (x$p!=a$p) stop("Argument dimensions do not match.\n")
  	if (x$p==1) {
		ITrgta <- ITrgt*as.matrix(a) 
		aITrgtU <- as.matrix(a) * ITrgtU
	}
  	else {
		ITrgta <- RightMult(ITrgt,as.matrix(a)) 
		aITrgtU <- as.matrix(a) %*% ITrgtU
	} 
  }
  else {
	if (is.matrix(a)) {
		if (x$p!=ncol(a)) stop("Argument dimensions do not match.\n")
        }
        else {
		if (x$p!=length(a))  stop("Argument dimensions do not match.\n")
		dim(a) <- c(x$p,1)
	}
  	if (x$p==1) {
		ITrgta <- matrix(ITrgt*drop(a),1,length(a)) 
		aITrgtU <- matrix(drop(a)*ITrgtU,length(a),1) 
	}
  	else {
		ITrgta <- RightMult(ITrgt,a) 
		aITrgtU <- t(a) %*% ITrgtU
	} 
   }
   if (x$q==1) {
	if (x$p==1) return( ITrgta - ITrgtU*aITrgtU/InnerM )
	return( ITrgta - ITrgtU%*%t(aITrgtU)/InnerM )
   }
   ITrgta - ITrgtU %*% solve(InnerM,t(aITrgtU)) 
}


