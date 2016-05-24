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


# = =============================================================== =
# =  Massive unit test to check all possible combinations of input  =
# =  parameters and make sure that the output matches the output    =
# =  from the serial version.                                       =
# = =============================================================== =

simplefunc <- function (data,indices){
    d <- data[indices]
    result <- mean(d)
    return(result)
}


test.simple <- function() {
  set.seed(88)
  a = boot(discoveries, simplefunc, 1001)
  set.seed(88)
	b = pboot(discoveries, simplefunc, 1001)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test simple default")
  
  a = boot(discoveries, simplefunc, 1001, simple=TRUE)
	b = pboot(discoveries, simplefunc, 1001, simple=TRUE)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
	b2 <- b[[2]]
	checkTrue(abs(mean(b2) - mean(discoveries)) < 1/sqrt(1001),"Test simple equals true")
  
  set.seed(88)
  a = boot(discoveries, simplefunc, 1001,  simple=FALSE)
  set.seed(88)
	b = pboot(discoveries, simplefunc, 1001, simple=FALSE)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test simple equals false")

}


