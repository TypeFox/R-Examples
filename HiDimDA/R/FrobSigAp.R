### FrobSigAp.R  (2011-06-13)
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

ForbSigap <- function(...) FrobSigAp(...)  # Wrapper, to ensure compability with version 0.1.0

FrobSigAp <- function(Sigma,q,nstarts=1,k0=NULL,penF=NULL,atol=1E-20,rtol=sqrt(.Machine$double.eps))
{
   p <- nrow(Sigma)
   if (p != ncol(Sigma)) stop("Matrix to be approximated (Sigma) is not square\n")

   if ( is.null(k0) || is.null(penF) ) {
  	if (is.null(k0)) k0 <- 0.01*min(diag(Sigma))
  	if (is.null(penF)) penF <- 100*max(diag(Sigma))
   }
   np <- q*(q+1)/2+(p-q)*q
   pu <- array(dim=np)

   SigSDec <- eigen(Sigma,symmetric=TRUE)
   cnt <- 0
   for (a in 1:q) {
	pu[(cnt+1):(cnt+p-a+1)] <- sqrt(SigSDec$values[a])*SigSDec$vectors[a:p,a]
      cnt <- cnt+p-a+1
   }
   psd <- rep(sqrt(SigSDec$values[1]/p),np)

   it <- 100000
   method <- "nlminb"

   optres <- RepLOptim(fr=f,gr=fgrad,inphess=fhess,parmean=pu,parsd=psd,nrep=nstarts,method=method,niter=it,nvar=p,q=q,
			Sigma=Sigma,k0=k0,penF=penF,atol=atol,rtol=rtol)

   if (optres$convergence > 0)  { 
	if (nstarts==1) {
		if (q>1) warning("Sigma approximation routine did not converge for a model with ",q," factors. Check the optimization results (available on the 'optres' field of the SigFq objects) to identify potentials problems, and/or try different starting points (using the 'nrep' argument) for the local optimization routine.")
		else warning("Sigma approximation routine did not converge for a model with 1 factor. Check the optimization results (available on the 'optres' field of the SigFq objects) to identify potentials problems, and/or try different starting points (using the 'nrep' argument) for the local optimization routine.")
	}
	else {
		if (q>1) warning("Sigma approximation routine did not always converge for a model with ",q," factors. Check the optimization results (available on the 'optres' field of the SigFq objects) to identify potentials problems.")
		else warning("Sigma approximation routine did not always converge for a model with 1 factor. Check the optimization results (available on the 'optres' field of the SigFq objects) to identify potentials problems.")
	}
   }

   B <- buildB(optres$par,p,q)
   if (p>1) D <- diag(Sigma) - apply(B,1,l2vnorm)
   else D <- Sigma[1,1] - apply(B,1,l2vnorm)
   SigFq(D,B,p,q,optres) # return(SigFq(D,B,p,q,optres))
}

buildB <- function(par,p,q)
{ 
   B <- matrix(0.,p,q)
   for (i in 1:q) B[i:p,i] <- par[((i-1)*(p-0.5*(i-2))+1):(i*(p-0.5*(i-1)))]
   B  # return(B)
}


