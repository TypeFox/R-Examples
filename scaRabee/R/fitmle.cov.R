
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

fitmle.cov <- function(problem=NULL,Fit=NULL){

# Preparing of fitted data (i.e., the observations)
  trts <- problem$data$trts
  ydata <- c()
  
  for (i in trts){
    ydata <- c(ydata,problem$data[[i]]$ana$DV)
  }
  
# Computation of the covariance matrix. The function also returns indices
  # to reorder x , in case the paramaters were not initially entered according
  # to P, L, IC, V order. It does the same to Fit.estimations for reporting
  # reasons.
  tmp <- get.cov.matrix(problem=problem,Fit=Fit)
  Fit <- c(Fit,list(cov=tmp$covmatrix,
                    orderedestimations=tmp$estimordered))
  rm(tmp)
  x <- transpose(Fit$orderedestimations$value)
  
  # Computes a series of derived statistics
    
  if (!any(Fit$cov=='singular')){
    # upper triangle of the correlation matrix
    Fit$cor <- Fit$cov/sqrt(transpose(diag(Fit$cov))%*%diag(Fit$cov))
    Fit$cor[lower.tri(Fit$cor)] <- 0
    
    # coefficients of variation
    Fit$cv  <- 100*sqrt(diag(Fit$cov))/x
    
    # interval used for confidence interval calculation
    delta   <- sqrt(diag(Fit$cov))*qt(0.975,length(ydata)-length(x))
    
    # confidence intervals
    Fit$ci  <- cbind(x-delta,x+delta)
      
  }
  
  # Aikaike information criterium
  Fit$AIC <- 2*(Fit$fval+length(x))
  
  # Computes statistics on secondary parameters, if any
  Fit$sec <- get.secondary(subproblem=problem,x=Fit$estimations)
  if (length(Fit$sec$estimates)!=0){
    if (!any(Fit$cov=='singular')){
      Fit$sec$cov <- transpose(Fit$sec$pder) %*% Fit$cov %*% Fit$sec$pder
      Fit$sec$cv  <- 100*sqrt(diag(Fit$sec$cov))/Fit$sec$estimates
      delta       <- sqrt(diag(Fit$sec$cov))*qt(0.975,length(ydata)-length(x))
      Fit$sec$ci  <- cbind(Fit$sec$estimates-delta,Fit$sec$estimates+delta)
    }
  } else {
    Fit$sec$cov <- NULL
    Fit$sec$cv  <- NULL
    Fit$sec$ci  <- NULL
  }  
  
  return(Fit)
  
}

