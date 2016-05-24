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

CL_P <-
function(x,y,s,mu,A,c,inds){
	n <- length(inds)
	k <- length(x)
	out <- matrix(0,k,n)
	for(i in 1:n){
		ki <- length(inds[[i]])
		w1bar <- x[inds[[i]]] + mu[i]*s[inds[[i]]]
		w2bar <- mu[i]*x[inds[[i]]] + s[inds[[i]]]
		wbar <- jordan_sqrt(jordan_square(w1bar) + jordan_square(w2bar) + 2*(mu[i]^2)*jordan_identity(ki))
		Lwinv <- solve(jordan_sym(wbar))
		out[inds[[i]],i] <- x[inds[[i]]] + s[inds[[i]]] - Lwinv %*% (jordan_sym(w1bar) %*% s[inds[[i]]] + jordan_sym(w2bar) %*% x[inds[[i]]] + 2*mu[i]*jordan_identity(ki))
	}
	return(out)
}

