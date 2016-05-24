### FrobSigAp1.R  (2011-06-13)
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

ForbSigap1 <- function(...) FrobSigAp1(...)  # Wrapper, to ensure compability with version 0.1.0

FrobSigAp1 <- function(SigmaSr,SigmaRank,q,nstarts=1,k0=NULL,penF=NULL,atol=1E-20,rtol=100*sqrt(.Machine$double.eps))
{
   variances <- function(j) return(SigmaSr[,j]%*%SigmaSr[,j])

   p <- ncol(SigmaSr)
   SigmaDiag <- sapply(1:p,variances)
   if (is.null(k0)) k0 <- 0.01*min(SigmaDiag)
   if (is.null(penF)) penF <- 100*max(SigmaDiag)
   np <- q*(q+1)/2+(p-q)*q
   pu <- array(dim=np)

   SigSrSVD <- rghtsngv(SigmaSr)
   cnt <- 0
   for (a in 1:q) {
	pu[(cnt+1):(cnt+p-a+1)] <- SigSrSVD$d[a]*SigSrSVD$v[a:p,a]
	cnt <- cnt+p-a+1
   }
  psd <- rep(SigSrSVD$d[1]/sqrt(p),np)

   it <- 100000
   method <- "nlminb"

   optres <- RepLOptim(fr=f1,gr=fgrad1,parmean=pu,parsd=psd,nrep=nstarts,method=method,niter=it,nvar=p,q=q,
			SigmaSr=SigmaSr,SigmaRank=SigmaRank,k0=k0,penF=penF,atol=atol,rtol=rtol)

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
   D <- SigmaDiag - apply(B,1,l2vnorm)
   SigFq(D,B,p,q,optres) # return(SigFq(D,B,p,q,optres))
}
