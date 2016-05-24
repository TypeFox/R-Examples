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



stypefm3 <- function(data, ind, them){
 y = data[ind]
return(c(mean(y), mean(them) ))
}


test.stypefm3 <- function() {
  set.seed(88)
  a = boot(discoveries, stypefm3, 1001, stype="f", m=1)
  set.seed(88)
	b = pboot(discoveries, stypefm3, 1001, stype="f", m=1)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test type = f and m is 1")

 set.seed(88)
  a = boot(discoveries, stypefm3, 1001, stype="f", m=2)
  set.seed(88)
	b = pboot(discoveries, stypefm3, 1001, stype="f", m=2)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test type = f and m is 2")

 set.seed(88)
  a = boot(discoveries, stypefm3, 1001, stype="f", m=4)
  set.seed(88)
	b = pboot(discoveries, stypefm3, 1001, stype="f", m=4)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test type = f and m is 4")

 set.seed(88)
  a = boot(discoveries, stypefm3, 1001, stype="f", m=16)
  set.seed(88)
	b = pboot(discoveries, stypefm3, 1001, stype="f", m=16)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test type = f and m is 1")

}


