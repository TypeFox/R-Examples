
#Copyright (c) 2009-2014 Sebastien Bihorel
#All rights reserved.
#
#This file is part of scaRabee.
#
#    scaRabee is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    scaRabee is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with scaRabee.  If not, see <http://www.gnu.org/licenses/>.
#

pder <- function(subproblem=NULL,x=NULL){
  
  # Evaluate the system at the point estimates
  tmp <- problem.eval(subproblem=subproblem,x=x)
    F2 <- apply(subproblem$data$data,1,
                function(x,...) {
                  tmp$f[x[2]+1,which(tmp$f[1,]==x[1])]},
                tmp)
    W2 <- apply(subproblem$data$data,1,
                function(x,...) {
                  tmp$weight[x[2]+1,which(tmp$weight[1,]==x[1])]},
                tmp)
  rm(tmp)
  
  # Remove instances of F2 and W2 corresponding to NaN values in ydata
  ydata <- subproblem$data$data$DV
  
  # Compute the partial derivatives
  mpder <- matrix(NA,nrow=length(x),ncol=length(ydata))
  wpder <- matrix(NA,nrow=length(x),ncol=length(ydata))
  
  for (i in 1:length(x)){
    h <- 1e-7
    x2 <- x
    x2[i] <- x2[i]+h
    
    # Evaluate the system at the point estimates + h
    tmp <- problem.eval(subproblem=subproblem,x=x2)
      F3 <- apply(subproblem$data$data,1,
                  function(x,...) {
                    tmp$f[x[2]+1,which(tmp$f[1,]==x[1])]},
                  tmp)
      W3 <- apply(subproblem$data$data,1,
                  function(x,...) {
                    tmp$weight[x[2]+1,which(tmp$weight[1,]==x[1])]},
                  tmp)
    rm(tmp)        
    
    # Calculate derivatives
    mpder[i,] <- (F3-F2)/h
    wpder[i,] <- (W3-W2)/h
    
  }
  
  varargout <- list(mpder=mpder,wpder=wpder)
  
  return(varargout)
  
}
