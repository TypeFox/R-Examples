# file blockmatrix3.R
#
# This file contains a script building a 2 x 2 block-matrix from 4 3 x 3 real matrices. 
# This script makes use of 'blockmatrix' function
#
# author: Emanuele Cordano on 22-02-2012

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

###############################################################################
rm(list=ls())
library(blockmatrix) 

A <- array(rnorm(9,mean=1),c(3,3))
B <- 0 #array(rnorm(9,mean=2),c(3,3))
C <- 0
D <- array(rnorm(9,mean=4),c(3,3))
F <- array(rnorm(9,mean=10),c(3,3))

M <- blockmatrix(names=c("A","0","D","0"),A=A,D=D,dim=c(2,2))
E <- blockmatrix(names=c("0","F","D","0"),F=F,D=D,dim=c(2,2))

R <- M+E
S <- solve(R)
P <- blockmatmult(R,E)
