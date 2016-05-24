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

myfunction <- function (data,indices){
    d <- data[indices]
    result <- mean(d)
    return(result)
}


test.standard <- function() {
  set.seed(88)
  a = boot(discoveries, myfunction, 100)
  set.seed(88)
	b = pboot(discoveries, myfunction, 100)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test 0")
}

test.standard_df <- function() {
  set.seed(88)
  a = boot(trees[,1], myfunction, 340)
  set.seed(88)
	b = pboot(trees[,1], myfunction, 340)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test 0")
}

test.standard_label <- function() {
  set.seed(88)
  a = boot(trees$Girth, myfunction, 340)
  set.seed(88)
	b = pboot(trees$Girth, myfunction, 340)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test data with labels")
}

test.standard_sample <- function() {
  set.seed(88)
  a = boot(c(9,4,563,2,3,4,66,53.4), myfunction, 340)
  set.seed(88)
	b = pboot(c(9,4,563,2,3,4,66,53.4), myfunction, 340)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Test data is expression")
}

