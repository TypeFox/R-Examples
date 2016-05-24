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

FB_smoother <-
function(x, s, mu, inds){
	#This function has been generalized to allow for n>1
#	phi <- numeric(length(unlist(inds)))
	k <- length(x)
	phi <- numeric(k)
	for(i in 1:length(inds)){
		phi[inds[[i]]] <- (1 + mu[i])*(x[inds[[i]]] + s[inds[[i]]]) - jordan_sqrt(jordan_square(x[inds[[i]]] + mu[i]*s[inds[[i]]]) + jordan_square(mu[i]*x[inds[[i]]] + s[inds[[i]]]) + 2*(mu[i]^2)*jordan_identity(length(inds[[i]])))
	}
	return(phi)

}

