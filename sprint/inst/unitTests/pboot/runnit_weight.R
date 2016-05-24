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



stypewfun <- function(d, w) sum(d$x * w)/sum(d$u * w)

stypewmfun <- function(d, w, mym){
 return(c( sum(d$x * w)/sum(d$u * w), sum(mym)))
}


test.stypew <- function() {
  set.seed(88)
  a = boot(city, stypewfun, 4638, stype="w")
  set.seed(88)
	b = pboot(city, stypewfun, 4638, stype="w")
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test type = w ")

}

test.stypewm <- function() {
  set.seed(88)
  a = boot(city, stypewmfun, 1638, stype="w", m=1)
  set.seed(88)
	b = pboot(city, stypewmfun, 1638, stype="w", m=1)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test type = w m=1")

  set.seed(88)
  a = boot(city, stypewmfun, 638, stype="w", m=2)
  set.seed(88)
	b = pboot(city, stypewmfun, 638, stype="w", m=2)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test type = w m=2")

  set.seed(88)
  a = boot(city, stypewmfun, 638, stype="w", m=8)
  set.seed(88)
	b = pboot(city, stypewmfun, 638, stype="w", m=8)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test type = w m=8")

  set.seed(88)
  a = boot(city, stypewmfun, 468, stype="w", m=16)
  set.seed(88)
	b = pboot(city, stypewmfun, 468, stype="w", m=16)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test type = w m=16")

}


