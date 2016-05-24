### CInterface.R  (2011-04-20)
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

f <- function(par,nvar,q,k0,penF,Sigma)
{
	.Call("RFnDist",
			as.double(par),
			as.integer(nvar),
			as.integer(q),
			as.double(Sigma[row(Sigma)>=col(Sigma)]),
			as.double(k0),
			as.double(penF) 
		)
}

fgrad <- function(par,nvar,q,k0,penF,Sigma)
{
	.Call("Rfgrad",
		as.double(par),
		as.integer(nvar),
		as.integer(q),
		as.double(Sigma[row(Sigma)>=col(Sigma)]),
		as.double(k0),
		as.double(penF)
	)
}

fhess <- function(par,nvar,q,k0,penF,Sigma)
{
	.Call("Rfhess",
		as.double(par),			
		as.integer(nvar),
		as.integer(q),
		as.double(Sigma[row(Sigma)>=col(Sigma)]),
		as.double(k0),
		as.double(penF)
	)
}

f1 <- function(par,nvar,q,k0,penF,SigmaSr,SigmaRank,lbound=NULL,ubound=NULL)
{
	if (any(!is.finite(par))) return(Inf)
	.Call("RFnDist1",
			as.double(par),
			as.integer(nvar),
			as.integer(q),
			as.double(SigmaSr),
			as.integer(SigmaRank),
			as.double(k0),
			as.double(penF) 
		)
} 

fgrad1 <- function(par,nvar,q,k0,penF,SigmaSr,SigmaRank,lbound=NULL,ubound=NULL)
{
	.Call("Rfgrad1",
		as.double(par),
		as.integer(nvar),
		as.integer(q),
		as.double(SigmaSr),
		as.integer(SigmaRank),
		as.double(k0),
		as.double(penF)
	)
}

fhess1 <- function(par,nvar,q,k0,penF,SigmaSr,SigmaRank)
{
	.Call("Rfhess1",
		as.double(par),			
		as.integer(nvar),
		as.integer(q),
		as.double(SigmaSr),
		as.integer(SigmaRank),
		as.double(k0),
		as.double(penF)
	)
}


