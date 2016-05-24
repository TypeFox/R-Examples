# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2012-2013 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
#                    Jianchang Lin and Stuart Lipsitz.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#######################################################################
## density function for the TBS
## public version checks if arguments are valid
## dist is a list, we assume the four functions of the dist are
## given as $d(x,xi) $p(x,xi) $r(x,xi) $q(x,xi) plus a test function to check if
## the list of parameters are in compliance with the functions. If dist is
## a string, then dist.error is called to try to transform it into the necessary
## functions
dtbs <- function(time,lambda=1,xi=1,beta=1,x=NULL,dist=dist.error("norm")) {
  if(is.character(dist)) dist=dist.error(dist)
  aux <- .test.tbs(lambda,xi,beta,x,time=time,type="d")
  return(.dtbs(time,lambda,xi,aux$beta,aux$x,dist))
}
## this private version does not check the arguments, but assumes they are ok
## calling the private version is faster than the public because of that
.dtbs <- function(time,lambda=1,xi=1,beta=1,x=NULL,dist="norm") {
  out <- ((1/time)*(abs(log(time))^(lambda-1))*
          dist$d(c(.g.lambda(log(time),lambda)-.g.lambda(c(x%*%beta),lambda)),xi))    
  return(c(out))
}
