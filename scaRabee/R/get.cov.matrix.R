
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

get.cov.matrix <- function(problem=NULL,Fit=NULL){
  
  # Initialization of the matrix
  M <- matrix(0,nrow=length(Fit$estimations),ncol=length(Fit$estimations))
  
  # Reorder param list
  ordered <- order.parms.list(x=problem$init)
  
  # Filter the param and ordered data structures to get only the estimated
  # parameters
  estparam <- problem$init[which(problem$init$isfix==0),]
  fixparam <- problem$init[which(problem$init$isfix==1),]
  estorder <- ordered[which(ordered$isfix==0),]
  fixorder <- ordered[which(ordered$isfix==1),]
  
  # Calculates the number of estimated model and variance parameters
    # p = nb of model parametres
  p <- length(get.parms.data(x=estparam,which='type',type='P')) +
       length(get.parms.data(x=estparam,which='type',type='L')) +
       length(get.parms.data(x=estparam,which='type',type='IC'))
    # q = nb of variance parametres
  q <- length(get.parms.data(x=estparam,which='type',type='V')) 
  
  # Determines in estparam the corresponding parameter indices from estorder
  indices <- match(estorder$names,estparam$names)
  
  # Computes model predictions and partial derivatives using original
  # parameter order (x and param)
  w <- c() ; mpder <- c() ; wpder <- c()
  
  trts <- problem$data$trts
  
  for (i in trts){
    # Creates subproblem
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
    
    # Obtains the model prediction and weighting based on final estimates
    tmp <- problem.eval(subproblem=subproblem,x=Fit$estimations)
      subf <- apply(subproblem$data$data,1,
                    function(x,...) {
                      tmp$f[x[2]+1,which(tmp$f[1,]==x[1])]},
                    tmp)
      subw <- apply(subproblem$data$data,1,
                    function(x,...) {
                      tmp$weight[x[2]+1,which(tmp$weight[1,]==x[1])]},
                    tmp)
    rm(tmp)
    
    # Obtains the matrix or partial derivatives with respect to final estimates
    tmp <- pder(subproblem=subproblem,x=Fit$estimations)
      submpder <- tmp$mpder
      subwpder <- tmp$wpder
    rm(tmp)
    
    # Reorders submpder and subwpder to compute and output the matrix properly (P,L,IC,V)
    submpder <- submpder[indices,]
    subwpder <- subwpder[indices,]
    
    # Appends w, submpder, subwpder
    w <- c(w,subw)
    mpder <- cbind(mpder,submpder)
    wpder <- cbind(wpder,subwpder)
  }
  
  # Computes the matrix to be inversed to get the covariance matrix
  M <- matrix(NA,nrow=p+q,ncol=p+q)
  
  if (p>0){
    for (j in 1:p){
      # covariance of the model parameters
      for (k in 1:p){
        M[j,k] <- 0.5*sum((1/w^2)*wpder[j,]*wpder[k,])+sum((1/w)*mpder[j,]*mpder[k,])
      }
      if (q>0){
        # covariance of the model parameters with the variance parameters
        for (k in (p+1):(p+q)){
          M[j,k] <- 0.5*sum((1/w^2)*wpder[j,]*wpder[k,])
        }
      }
    }
  }
  
  if (q>0){
    for (j in (p+1):(p+q)){
      # covariance of variance parameters
      for (k in (p+1):(p+q)){
        M[j,k] <- 0.5*sum((1/w^2)*wpder[j,]*wpder[k,])
      }
      if (p>0){
        # covariance of the variance parameters with the model parameters
        for (k in 1:p){
          M[j,k] <- M[k,j]
        }
      }
    }
  }
  
  # Computes the covariance matrix
  if (det(M)==0){
    covmatrix <- 'singular'
  } else {
    covmatrix <- solve(qr(M,LAPACK=TRUE))
  }
  
  # Gets estimated parameter reordered values, names, type
  estimordered <- data.frame(names=estparam$names[indices],
                             value=Fit$estimations[indices],
                             type =estparam$type[indices],
                             stringsAsFactors=FALSE)
  
  varargout <- list(covmatrix=covmatrix,estimordered=estimordered)
  
  return(varargout)
  
}

