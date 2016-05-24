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

CL_H <-
function(x, s, mu, A, b, inds){
#	S <- b - A %*% x
#	FB <- matrix(FB_smoother(x, s, mu, inds),,1)
#	mu <- mu
#	save(S,FB,mu,file="~/sandbox/SFBmutest.rda")
	n <- length(inds)
	k <- ncol(A)
	m <- nrow(A)
	out <- numeric(m+k+n)
	out[1:m] <- as.vector(b - A %*% x)
	out[(m+1):(m+k)] <- FB_smoother(x, s, mu, inds)
	out[(m+k+1):(m+k+n)] <- mu

#out <- c(as.matrix(b - A %*% x), FB_smoother(x, s, mu, inds), mu)

	return(out)
}

