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

convert_l <-
function(kvec, type){
	newkvec <- integer()
	newtype <- character()
	for(i in 1:length(type)){
		if (type[i] == 'l'){
			newkvec <- c(newkvec,rep(1,kvec[i]))
			newtype <- c(newtype,rep('q',kvec[i]))
		}else if(type[i] == 'q'){
			newkvec <- c(newkvec,kvec[i])
			newtype <- c(newtype,type[i])
		}
		
	}
	return(list(kvec=newkvec,type=newtype))
}

