
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

residual.report <- function(problem=NULL,Fit=NULL,files=NULL){
  
  # Definition of some local variables by transfer from Fit and states variables
    # Optimized model and variance parameters
  x <- Fit$estimations
  
  # Copies param in newparam and replaces the new estimates in newparam
  newparam <- problem$init
  npar     <- length(newparam$names)
  estindex <- 1
  for (i in 1:npar){
    if (newparam$isfix[i]==0){
      newparam$value[i] <- x[estindex]
      estindex <- estindex + 1
    }
  }
  
  # Consider each dose level
  trts <- problem$data$trts
  
  for (i in trts){
    # Creates subproblem
    subproblem <- problem[c('code','method','init','debugmode','modfun',
                            'solver.options')]
    subproblem$data$xdata <- sort(unique(problem$data[[i]]$ana$TIME))
    ana.data <- problem$data[[i]]$ana
    subproblem$bolus <- problem$data[[i]]$bolus
    subproblem$infusion <- problem$data[[i]]$infusion
    
    if (size(problem$data[[i]]$cov,1)!=0){
      subproblem$cov <- problem$data[[i]]$cov
    } else {
      subproblem$cov <- list(NULL)
    }
    
    # Get the model predictions and corresponding weights
    pred <- problem.eval(subproblem,x)
    ana.data$IPRED <- apply(ana.data,1,
                            function(x,...) {
                              pred$f[x[2]+1,which(pred$f[1,]==x[1])]},
                            pred)
    ana.data$W <- apply(ana.data,1,
                        function(x,...) {
                          pred$weight[x[2]+1,which(pred$weight[1,]==x[1])]},
                        pred)
    
    # Calculates the residuals
    ana.data$IRES <- ana.data$DV - ana.data$IPRED
    ana.data$WRES <- ana.data$IRES/ana.data$W
    
    # Add variables to ana.data
    ana.data$ID <- problem$data$id
    ana.data$TRT <- i
    
    # Print the residuals to file
    tmpform <- paste(c('%d,%s,%d,%d,',rep('%0.5g,',5),'%0.5g'),collapse='')
    for (row in 1:size(ana.data,1)){
      tmp <- c(ana.data$ID[row],ana.data$TRT[row],row,ana.data$DVID[row],
               ana.data$TIME[row],ana.data$DV[row],ana.data$IPRED[row],
               ana.data$IRES[row],ana.data$W[row],ana.data$WRES[row])
      write(do.call(sprintf,c(list(tmpform),as.list(tmp))),
            file=files$pred,append=TRUE,sep='\n')
    }
  }
}

