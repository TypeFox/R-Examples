
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

scarabee.check.model <- function(problem=NULL,files=NULL){
  
  # Check inputs
  if (is.null(problem) | is.null(files)){
    stop(paste('one or more input argument of fitmle is null. Please, ',
               'check your code.',sep=''))
  }
  
  # Check model
  trts <- problem$data$trts
  
  for (i in trts){
    # Create subproblem
    subproblem             <- problem[c('code','method','init','debugmode',
                                        'modfun','solver.options')]
    subproblem$data$xdata  <- sort(unique(problem$data[[i]]$ana$TIME))
    subproblem$data$data  <- problem$data[[i]]$ana
    subproblem$bolus <- problem$data[[i]]$bolus
    subproblem$infusion <- problem$data[[i]]$infusion
    
    if (size(problem$data[[i]]$cov,1)!=0){
      subproblem$cov <- problem$data[[i]]$cov
    } else {
      subproblem$cov <- list(NULL)
    }
    
    # vector of estimated parameters
    x <- problem$init[which(problem$init$isfix==0),'value']
    
    # Get the model predictions and corresponding weights
    pred <- problem.eval(subproblem=subproblem,x=x,check=TRUE)
  }
}
