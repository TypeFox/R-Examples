##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################


tnfunction <- function (data, indices){
    girth <- data[indices, 1]
    height <- data[indices, 2]
    volume <- data[indices, 3]
    return(c(mean(girth), mean(height), mean(volume)))
} 

test.tn <- function() {
  
  set.seed(88)
  a = boot(trees, tnfunction, 1000)
  set.seed(88)
	b = pboot(trees, tnfunction, 1000)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"test tn")
}


