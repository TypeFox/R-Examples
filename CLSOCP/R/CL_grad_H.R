######################################################################
#Copyright Jason Rudy & Faramarz Valafar 2009-2010

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
######################################################################

CL_grad_H <-
function(x, y, s, mu, A, c, inds){
	m <- dim(A)[1]
	k <- dim(A)[2]
	n <- length(inds)

	if (get("use_sparse",pos=parent.frame(1))){
		out <- Matrix(0, m+k+n, m+k+n)
#	print("Using sparse Jacobian")
	}else{
		out <- matrix(0, m+k+n, m+k+n)
	}

	
	out[1:m,1:k] <- -1*A
	

	N <- matrix(0,k,k)
	P <- matrix(0,k,n)
	M <- matrix(0,k,k)
	for(i in 1:n){
		ki <- length(inds[[i]])
		w1bar <- x[inds[[i]]] + mu[i]*s[inds[[i]]]
		w2bar <- mu[i]*x[inds[[i]]] + s[inds[[i]]]
		wbar <- jordan_sqrt(jordan_square(w1bar) + jordan_square(w2bar) + 2*(mu[i]^2)*jordan_identity(ki))
		Lwinv <- solve(jordan_sym(wbar))
		M[inds[[i]],inds[[i]]] <- (1+mu[i])*diag(ki) - Lwinv %*% (jordan_sym(w1bar) + mu[i]*jordan_sym(w2bar))
		N[inds[[i]],inds[[i]]] <- (1+mu[i])*diag(ki) - Lwinv %*% (mu[i]*jordan_sym(w1bar) + jordan_sym(w2bar))
		P[inds[[i]],i] <- x[inds[[i]]] + s[inds[[i]]] - Lwinv %*% (jordan_sym(w1bar) %*% s[inds[[i]]] + jordan_sym(w2bar) %*% x[inds[[i]]] + 2*mu[i]*jordan_identity(ki))
		
	}

	
	
	
	out[(m+1):(m+k),1:k] <- M #CL_M(x,y,s,mu,A,c,inds)

	out[(m+1):(m+k),(k+1):(k+m)] <- -1*N %*% t(A)
	out[(m+1):(m+k),(k+m+1):(k+m+n)] <- P #CL_P(x,y,s,mu,A,c,inds)
	out[(m+k+1):(m+k+n),(m+k+1):(m+k+n)] <- diag(n)
	

	return(out)
}

