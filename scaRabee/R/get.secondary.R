
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

get.secondary <- function(subproblem=NULL,x=NULL){
  
  # Computes secondary parameter estimates
  secparam <- compute.secondary(subproblem=subproblem,x=x)
  
  # Computes partial derivatives of secondary parameters
  secparam$pder <- matrix(NA,nrow=length(x),ncol=length(secparam$estimates))
  if (length(secparam$estimates)!=0){
    for (i in 1:length(x)){
      h <- 1e-7 # RTOL ATOL from ode solver
      
      x2 <- x
      x2[i] <- x2[i]+h
      
      tmp <- compute.secondary(subproblem=subproblem,x=x2)
        sec2 <- tmp$estimates
      rm(tmp)
      secparam$pder[i,] <- transpose((sec2-secparam$estimates)/h)
    }
    
    # Reorders param data structure
    ordered <- order.parms.list(x=subproblem$init)
    
    # Filters the param and ordered data structures to get only the estimated
    # parameters
    estparam <- subproblem$init[which(subproblem$init$isfix==0),]
    fixparam <- subproblem$init[which(subproblem$init$isfix==1),]
    estorder <- ordered[which(ordered$isfix==0),]
    fixorder <- ordered[which(ordered$isfix==1),]
    
    # Calculates the number of estimated model and variance parameters
      # p = nb of model parametres
    p <- length(get.parms.data(x=estparam,which='type',type='P'))+ 
         length(get.parms.data(x=estparam,which='type',type='L'))+ 
         length(get.parms.data(x=estparam,which='type',type='IC'))
      # q = nb of variance parametres
    q <- length(get.parms.data(x=estparam,which='type',type='V'))
    
    # Determines in estparam the corresponding parameter indices from estorder
    indices <- c()
    for (i in 1:(p+q)){
      for (j in 1:(p+q)){
         if (estparam$names[j]==estorder$names[i]){
           indices <- c(indices,j)
         }
      }
    }
    
    # Reorders secparam.pder
    secparam$pder <- secparam$pder[indices,,drop=FALSE]
  }
  
  return (secparam)
  
}

