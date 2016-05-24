###########################################################################
# Copyright 2009 Nobody                                                   #
#                                                                         #
# This file is part of hergm.                                             #
#                                                                         # 
#    hergm is free software: you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         # 
#    hergm is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
#    GNU General Public License for more details.                         #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with hergm.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                         # 
###########################################################################

hergm.permutation.wrapper <- function(number) 
{
  number_permutations <- factorial(number)
  permutations <- vector(length=number_permutations*number)
  for (i in 1:number) permutations[i] = i
  output <- .C("Permutations",
                     as.integer(number),
                     as.integer(number_permutations),
                     permutations = as.integer(permutations),
                     PACKAGE="hergm")
  permutations <- matrix(output$permutations, nrow = number_permutations, ncol = number, byrow = TRUE)
  permutations
}

