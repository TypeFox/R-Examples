# 		RevNet: tools for reverse engineering of networks, part of the 'MetStaT' package  
#		Copyright (C) 2012 Diana Hendrickx and Tim Dorscheidt
#		
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation, either version 3 of the License, or
#		(at your option) any later version.
#		
#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.
#		
#		You should have received a copy of the GNU General Public License
#		along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# 		Email: g.zwanenburg@uva.nl ('MetStaT' contact person) or tdorscheidt@gmail.com
#		Article with further details:
#		Reverse engineering of metabolic networks, a critical assessment
#		Diana M. Hendrickx, Margriet M. W. B. Hendriks, Paul H. C. Eilers, Age K. Smilde and Huub C. J. Hoefsloot 
#		Mol. BioSyst, Volume 7:2 (2011) pages 511-520
###############################################################################

### Three different approaches for reverse engineering of networks are included in this file, as well as some more general methods at the end.

## Reverse Engineering of Networks using the 'Jacobian' approach
RevNet.JacobianMethod <- function(data, delta.t, steady.state.concentrations, lamba.penalty.par, kappa.penalty.par, jacobian.threshold) {
	# Internal method for the Jacobian method
	RevNet.CalcSSQ <- function(matrix) {
		sum(matrix^2)
	}
	
	# Internal method for the Jacobian method	
	RevNet.FourthOrder <- function(X, deltat) {
		dimX <- dim(X)
		Y <- array(dim=c(dimX[1]-4,dimX[2],dimX[3]))
		Y[,,] <- 
				( -X[5:dimX[1],,,drop=FALSE] 
					+ 8*X[4:(dimX[1]-1),,,drop=FALSE] 
					- 8*X[2:(dimX[1]-3),,,drop=FALSE] 
					+ X[1:(dimX[1]-4),,,drop=FALSE] ) / (12 * deltat)
		Y
	}
	
	# Internal method for the Jacobian method	
	RevNet.FirstEstimate <- function(X, Y, s) {
		dimX <- dim(X)
		n <- dimX[1]
		m <- dimX[2]
		p <- dimX[3]
		Xnew <- X[3:(dimX[1]-2),,,drop=FALSE]
		Xblock <- array(dim=c((n-4)*p,m))
		Yblock <- array(dim=c((n-4)*p,m))
		for (i in 1:p) {
			Xblock[( (1+(i-1)*(n-4)) : (i*(n-4))) ,] = Xnew[,,i] - matrix(s,nrow=n-4,ncol=length(s),byrow=T)
			Yblock[( (1+(i-1)*(n-4)) : (i*(n-4))) ,] = Y[,,i,drop=FALSE]
		}
		Xfirst <- array(dim=c((n-4)*p*m,m*m))
		yfirst <- array(dim=c((n-4)*p*m,1))  
		for (j in 1:m) {
			Xfirst[( (1+(j-1)*(n-4)*p) : (j*(n-4)*p) ) , ( (1+(j-1)*m) : (j*m) )] = Xblock
			yfirst[( (1+(j-1)*(n-4)*p) : (j*(n-4)*p) ) , 1] = Yblock[,j,drop=FALSE];
		}
		Xfirst[is.na(Xfirst)] <- 0
		yfirst[is.na(yfirst)] <- 0
		first_train <- data.frame(Xfirst = Xfirst, yfirst = yfirst)
		R <- coef(lm.ridge(yfirst~., first_train, lambda = 0.001))
		jfirst <- R[2:(m^2+1) , drop=FALSE];
		JACfirst <- t(matrix(jfirst,nrow=m,ncol=m))
		list(JACfirst = JACfirst, jfirst = jfirst, Xfirst = Xfirst, yfirst = yfirst)
	}
	
	# Internal method for the Jacobian method	
	RevNet.LOL1reg <- function(jfirst, Xfirst, yfirst, lambda, kappa) {
		jnew <- jfirst;
		XX <- t(Xfirst) %*% Xfirst
		XY <- t(Xfirst) %*% yfirst
		factor <- lambda^2
		m <- length(jnew)^0.5
		a <- which(diag(m)==0)
		oonneess <- rep(1,length(jfirst[a]))
		v <- rep(0,length(jnew))
		v[a] = oonneess/(1e-6+abs(jnew[a])+kappa*jnew[a]^2)
		L <- diag(v) * factor
		jold = jnew
		jnew = MetStaT.mldivide((XX + L),XY)
		while ((RevNet.CalcSSQ(jold/jold - jnew/jold))>0.00001) {
			v = rep(0,length(jnew));
			v[a] = oonneess/(1e-6+abs(jnew[a])+kappa*jnew[a]^2)
			L = factor*diag(v)
			jold = jnew
			jnew = MetStaT.mldivide((XX + L),XY)
		}
		JAC = t(matrix(jnew,nrow=m,ncol=m))
		JAC
	}
	
	# Internal method for the Jacobian method	
	RevNet.VertexEdge <- function(JAC, threshold) {
		N <- matrix(0,nrow=dim(JAC)[1],ncol=dim(JAC)[2])
		N[abs(JAC)>=threshold] <- 1
		N
	}
	
	Y <- RevNet.FourthOrder(data, delta.t)
	res <- RevNet.FirstEstimate(data, Y, steady.state.concentrations)
	jac <- RevNet.LOL1reg(res$jfirst, res$Xfirst, res$yfirst, lamba.penalty.par, kappa.penalty.par)
	N <- RevNet.VertexEdge(jac, jacobian.threshold)
	N
}

## Reverse Engineering of Networks using the 'time lag' approach
RevNet.TimeLaggedMethod <- function(data, max.time.lag, threshold) {
	# Internal method for the Time Lagged method	
	RevNet.CcfMax <- function(X, maxlag, plotCor = FALSE) {
		dimX <- dim(X)
		n <- dimX[1]
		m <- dimX[2]
		
		C <- matrix(nrow = m, ncol = m)
		
		for (i in 1:m) {
			for (j in 1:m) {
				crcor <- ccf(X[,i], X[,j], lag.max = maxlag, type = c("correlation"), plot = plotCor)$acf
				C[i,j] <- max(abs(crcor))
			}
		}
		C
	}
	# Internal method for the Time Lagged method
	RevNet.Connection <- function(C, percent) {
		m <- dim(C)[1]
		C.lower.tri <- C[lower.tri(C, diag = FALSE)] # returns part below the diagonal
		C.filtered <- C.lower.tri[C.lower.tri!=0]
		p <- quantile(C.filtered, percent);

		N <- matrix(NA, nrow=m, ncol=m)
		for (i in 1:m) {
			for (j in 1:m) {
				if (i==j) {
					N[i,j]=1;
				} else if (C[i,j]<p) {
					N[i,j]=0;
				} else {
					N[i,j]=1;
				}
			}
		}
		N
	}
	
	res <- RevNet.CcfMax(data, max.time.lag, plotCor=FALSE)
	N <- RevNet.Connection(res, threshold)
	N
}

## Reverse Engineering of Networks using the 'zero slopes' approach
RevNet.ZeroSlopesMethod <- function(X, deltat, threshold) {
	n <- dim(X)[1]
	m <- dim(X)[2]
	p <- dim(X)[3]
	
	N <- matrix(NA, nrow=m, ncol=p)
	
	for (i in 1:p) {
		for (j in 1:m) {
			if (i==j) {
				N[j,i]=-1
			} else if ( ((abs(X[2,j,i]-X[1,j,i]))/deltat)<threshold ) {
				N[j,i]=0
			} else {
				N[j,i]=1
			}
		}
	}
	N
}
