######################################################################
#
# BayesH.R 
#
# copyright (c) 2016-2016, Renato Rodrigues Silva
# last modified Mar, 2016
# first written Jan, 2016
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Part of the BayesH package
# Contains:  bayesModel.fit     
#
######################################################################



########################################################################
#Function
########################################################################

bayesModel.fit = function(X, y, nu0, s0,  niter=2000, burnin=500, type ="bayesH")
{
    #Missing arguments
    if(missing(y)){
       stop("You must to provide a response variable")
    }
    if(missing(X)){
       stop("You must to provide an incidence matrix")
    }      
    #Checking errors
    if(!is.matrix(X) ){
        stop("X must be a matrix")
    } 
    if(any(is.na(X))){
        stop("NA's are not allowed in the incidence matrix")
    }
    if(!is.vector(y)  ){
        stop("y must be a vector")
    } 
    if(!is.vector(niter) || !is.vector(burnin) ){
        stop("These objects must be a vector")
    } 
    if(!( (type == "ridge")  || (type == "bayesH") ) ) {
       stop("Wrong type!")
    } 
    
    

    if(type == "ridge"){
      if(missing(nu0) ){
         stop("You must to provide the vector mu0")
      } 
      if(missing(s0)){
         s0 = get.scale.bayesH(nu0 = nu0, R2 = 0.5, X = X, y = y, type ="ridge")
      } 
      
      if( ( length(nu0) != 2) ||  (length(s0) != 2)  || !is.numeric(nu0)  || !is.numeric(s0) ){
          stop("The length of the vectors nu0 and s0 must be one and these vector must be a numeric.")
      }      
      if(any(c(nu0,s0) < 0)){
         stop("The values of the nu0 and s0 must be greater than zero")
      }
           
       
    } else if(type == "bayesH"){
      if(missing(nu0) ){
         stop("You must to provide the vector nu0")
      }         
      if((length(nu0) != 3) || (length(s0) != 3)  || !is.numeric(nu0) || !is.numeric(s0)){
         stop("The length of the vectors nu0 and s0 must be three and these vector must be a numeric.")
      }

      if(any(c(nu0,s0) < 0)){
         stop("All elements of the nu0 and s0 must be greater than zero")
      }

      if(missing(s0)){
         s0 = get.scale.bayesH(nu0 = nu0, R2 = 0.5, X = X, y = y, type = "bayesH")
      }           
      
    } else {
       stop("Error")
    }
    if((burnin %% 1 != 0) || (niter %% 1 != 0) ){
       stop("These numbers must be integer") 
    }        
    if(nrow(X) <= 2 || ncol(X) <= 2){
       stop("Data not available to fit the model")
    }     
     
    group = as.integer(is.na(y))
    ymis = y  
    if(any(is.na(y))) {       
       ymis[group == 1] = mean(y, na.rm = TRUE)
    }
    parBeta = matrix(rep(0, ncol(X)), ncol=1)  
    mu = mean(y, na.rm=TRUE)
    sig = var(y, na.rm=TRUE) 
    if(type == "ridge"){    
         bayesfit = .Call("bayesRidge", X, matrix(ymis, ncol=1), mu, parBeta,
                     as.double(sig),
                     as.integer(group),                                               
                     as.double(nu0),
                     as.double(s0),                      
                     as.integer(niter),
                     as.integer(burnin))  
    } else {          
         bayesfit = .Call("bayesH", X, matrix(ymis, ncol=1), mu, parBeta,
                     as.double(sig), 
                     as.integer(group),                                                  
                     as.double(nu0),
                     as.double(s0),                      
                     as.integer(niter),
                     as.integer(burnin))  
    }            
    sigmaGibbs = bayesfit[[1]]
    muGibbs = bayesfit[[2]]
    betaGibbs = bayesfit[[3]] 
    predGibbs = bayesfit[[4]] 
    #Output
    out = list( muGibbs = muGibbs, sigmaGibbs = sigmaGibbs,
                coef = betaGibbs, fitted.values = predGibbs ,y=y)                       
    class(out) = "BayesH"
    return(out)    
}





