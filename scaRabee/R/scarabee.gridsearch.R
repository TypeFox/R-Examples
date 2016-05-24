
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

scarabee.gridsearch <- function(problem=NULL,
                                npts=NULL,
                                alpha=NULL,
                                files=NULL) {
  
  # Check inputs
  if (is.null(problem) | is.null(npts) | is.null(files)){
    stop('one or more argument is null. Please, check your code.')
  }
  
  # Create objection function minimization routine
  grid.obj.fun <- function(x=NULL,...){
    
    # Get model predictions and calculates objective function (ML)
    ML <- 0
    
    for (id in problem$data$ids){
      
      # Create subproblem for id
      idproblem <- problem[c('code','method','init','debugmode','modfun',
                             'solver.options')]
      idproblem$data <- problem$data[[id]]
      
      # Compute ML for each id
      for (i in idproblem$data$trts){
        
        # Create subproblem
        subproblem <- idproblem[c('code','method','init','debugmode','modfun',
                                  'solver.options')]
        subproblem$data$xdata <- sort(unique(idproblem$data[[i]]$ana$TIME))
        subproblem$data$data <- idproblem$data[[i]]$ana
        subproblem$bolus <- idproblem$data[[i]]$bolus
        subproblem$infusion <- idproblem$data[[i]]$infusion
        
        if (size(idproblem$data[[i]]$cov,1)!=0){
          subproblem$cov <- idproblem$data[[i]]$cov
        } else {
          subproblem$cov <- list(NULL)
        }
        
        # Calculate ML
        if (!idproblem$debugmode){
          tmp <- try(
            {
            # Get the model predictions and corresponding weights
            pred <- problem.eval(subproblem,x,grid=TRUE)
            ydata <- subproblem$data$data$DV
            fpred <- apply(subproblem$data$data,1,
                           function(x,...) {
                             pred$f[x[2]+1,which(pred$f[1,]==x[1])]},
                           pred)
            weight <- apply(subproblem$data$data,1,
                            function(x,...) {
                              pred$weight[x[2]+1,which(pred$weight[1,]==x[1])]},
                            pred)
            # Calculate minus twice the log likelihood function
            ML <- ML+sum(0.5*((((fpred-ydata)^2)/(weight))+log(weight)+log(2*pi)))
            },
            silent=TRUE)
          
          if (class(tmp)=="try-error"){
            cat(paste('An error occured during the computation of the model',
                'predictions or the negative log likelihood.\n'))
            ML <- Inf
          } else {
            ML <- tmp
          }
            
        } else {
          # Get the model predictions and corresponding weights
          pred <- problem.eval(subproblem,x,grid=TRUE)
          ydata <- subproblem$data$data$DV
          fpred <- apply(subproblem$data$data,1,
                         function(x,...) {
                           pred$f[x[2]+1,which(pred$f[1,]==x[1])]},
                         pred)
          weight <- apply(subproblem$data$data,1,
                          function(x,...) {
                            pred$weight[x[2]+1,which(pred$weight[1,]==x[1])]},
                          pred)
          # Calculate minus twice the log likelihood function
          ML=ML+sum(0.5*((((fpred-ydata)^2)/(weight))+log(weight)+log(2*pi)))
        }
      }
      
      if (is.null(ML) | is.na(ML)){
        warning(paste('\nIn grid.obj.fun: likelihood function was NA or NULL; it was ',
                      'thus coerced to Inf.\nPlease, check your model',
                      ' and your parameter definition (especially lower and\n ',
                      'upper boundaries) for potential divisions by 0.',sep=''),
                call.=FALSE,
                immediate.=TRUE)
        ML <- Inf
      }
    }
    
    return(ML)
    
  }
  
  # Definition of parameter initials and bounds
  x0 <- problem$init$value
  xmin <- problem$init$lb
  xmax <- problem$init$ub
  
  if (!is.null(alpha)){
    index <- which(problem$init$isfix==0)
    xvar0 <- x0[index]
    if (length(alpha)>length(xvar0)){
      alpha <- alpha[1:length(xvar0)]
    } else {
      alpha <- rep(alpha, length.out=length(xvar0))
    }
    xmin[index] <- (xvar0/alpha)
    xmax[index] <- (xvar0*alpha)
  }
  
  index <- which(problem$init$isfix==1)
  xmin[index] <- xmax[index] <- x0[index]
  
  # initial guess for parameters to be fitted
  fgrid <- fmin.gridsearch(fun=grid.obj.fun,x0=x0,
                           xmin=xmin,xmax=xmax,
                           npts=npts,alpha=alpha)
  
  # rename fgrid columns
  names(fgrid)[1:length(x0)] <- problem$init$names 
  
  return(fgrid)
  
}

