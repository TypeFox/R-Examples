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

make_inds <-
function(kvec, type){
	tot <- sum(kvec)
	len <- length(kvec)
	cur <- 1
	block <- 1
	inds <- list()
	while(block <= len){
		inds[[block]] <- cur:(cur+kvec[block] - 1)
		cur <- cur + kvec[block]
		block <- block + 1
	}
	return(inds)
}

