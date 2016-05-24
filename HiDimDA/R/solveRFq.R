### SolveRFq.R  (2011-06-13)
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

solve.SigFq <- function(a,b=NULL,...)
{
  if (a$q==1) DiB <-  drop(a$B / a$D)
  else DiB <-  a$B / a$D

  if (a$p==1)  {
  	if (a$q==1) {
		InnerM <- 1. + a$B*DiB
		if (is.null(b)) {
			res <- list(D=1./a$D,B=matrix(DiB/sqrt(InnerM),a$p,a$q),p=a$p,q=a$q)
			class(res) <- "SigFqInv"
			return(res)
		}
		else {
			Dib <-  drop(b/a$D)
			return( Dib - (DiB*a$B/InnerM)*Dib )
		}
	}
  	else {
		InnerM <- diag(a$q) + drop(a$B)%o%drop(DiB)
		if (is.null(b)) {
			res <- list(D=1./a$D,B=t(matrix(solve(t(chol(InnerM)),t(DiB)),a$p,a$q)),p=a$p,q=a$q)
			class(res) <- "SigFqInv"
			return(res)
		}
		else {
			if (!is.matrix(b)) bcol  <- 1
			else bcol <- ncol(b)
			Dib <-  drop(b/a$D)
			return( Dib - DiB%*%matrix(solve(InnerM,t(a$B))%*%Dib,a$q,bcol) )
		}
	}
  }
  else  {
  	if (a$q==1) {
		InnerM <- 1. + drop(t(a$B)%*%DiB)
		if (is.null(b)) {
			res <- list(D=1./a$D,B=matrix(DiB/sqrt(InnerM),a$p,a$q),p=a$p,q=a$q)
			class(res) <- "SigFqInv"
			return(res)
		}
		else {
			Dib <-  b / a$D
			return( Dib - DiB*(t(a$B)%*%Dib)/InnerM )
		}
	}
  	else {
		InnerM <- diag(a$q) + t(a$B)%*%DiB
		if (is.null(b)) {
			res <- list(D=1./a$D,B=t(solve(t(chol(InnerM)),t(DiB))),p=a$p,q=a$q)
			class(res) <- "SigFqInv"
			return(res)
		}
		else {
			if (!is.matrix(b)) bcol  <- 1
			else bcol <- ncol(b)
			Dib <-  b / a$D
			return( Dib - DiB%*%matrix(solve(InnerM,t(a$B))%*%Dib,a$q,bcol) )
		}
  	}
  }
}

solve.SigFqInv <- function(a,b=NULL,...)
{
  if (a$q==1) DiB <-  drop(a$B / a$D)
  else DiB <-  a$B / a$D
  if (a$p==1)  {
  	if (a$q==1) {
		InnerM <- 1. - a$B*DiB
		if (is.null(b)) return(SigFq(D=1./a$D,B=matrix(DiB/sqrt(InnerM),a$p,a$q),p=a$p,q=a$q))
		else {
			Dib <-  drop(b/a$D)
			return( Dib + (DiB*a$B/InnerM)*Dib )
		}
	}
  	else {
		InnerM <- diag(a$q) - drop(a$B)%o%drop(DiB)
		if (is.null(b)) return(SigFq(D=1./a$D,B=t(matrix(solve(t(chol(InnerM)),t(DiB)),a$p,a$q)),p=a$p,q=a$q))
		else {
			if (!is.matrix(b)) bcol  <- 1
			else bcol <- ncol(b)
			Dib <-  drop(b/a$D)
			return( Dib + DiB%*%matrix(solve(InnerM,t(a$B))%*%Dib,a$q,bcol) )
		}
	}
  }
  else  {
  	if (a$q==1) {
		InnerM <- 1. - drop(t(a$B)%*%DiB)
		if (is.null(b)) return(SigFq(D=1./a$D,B=matrix(DiB/sqrt(InnerM),a$p,a$q),p=a$p,q=a$q))
		else {
			Dib <-  b / a$D
			return( Dib + DiB*(t(a$B)%*%Dib)/InnerM )
		}
	}
  	else {
		InnerM <- diag(a$q) - t(a$B)%*%DiB
		if (is.null(b)) return(SigFq(D=1./a$D,B=t(solve(t(chol(InnerM)),t(DiB))),p=a$p,q=a$q))
		else {
			if (!is.matrix(b)) bcol  <- 1
			else bcol <- ncol(b)
			Dib <-  b / a$D
			return( Dib + DiB%*%matrix(solve(InnerM,t(a$B))%*%Dib,a$q,bcol) )
		}
  	}
  }
}
