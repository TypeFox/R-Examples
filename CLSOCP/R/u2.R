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

u2 <-
function(x){
	nrm <- two_norm(x[-1])
	n <- length(x)
	if(nrm < 1e-12){
		c(0.5,0.5*rep(sqrt(1/(n-1)),n-1))
	}else{
		c(0.5, 0.5*(1/nrm)*x[-1])
	}
}

