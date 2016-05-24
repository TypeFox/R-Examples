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
# Contains:  get.scale.bayesH  
#
######################################################################



########################################################################
#Function
########################################################################


get.scale.bayesH = function(nu0, R2 = 0.5, X, y, type = "bayesH")
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
    if(missing(nu0) ){
       stop("You must to provide the value of nu0")
    } 
    if(any(c(nu0) < 0)){
       stop("All elements of the nu0  must be greater than zero")
    }
    if(R2 > 1 || R2 < 0 ){
      stop("r squared must be between zero and one") 
    }

    if(type=="ridge"){
      if(length(nu0) != 2){
         stop("The length of nu0 must be two")
      } 
 
      s0 = c(NA, NA)
      s0[2] = ( (nu0[2] + 2) / nu0[2] ) * (1 - R2) * var(y, na.rm = TRUE)        
      xx = apply(X, 2, function(x){ sum(x^2, na.rm = TRUE) })
      s0[1] = (1 / mean(xx) ) * R2 * var(y, na.rm = TRUE) * ( (nu0[1] + 2) / nu0[1])
      
    } else if(type == "bayesH"){
      if(length(nu0) != 3){
         stop("The length of nu0 must be three")
      } 
 
      s0 = c(NA, NA, NA)
      s0[3] = ( (nu0[3] + 2) / nu0[3] ) * (1 - R2) * var(y, na.rm = TRUE)        
      xx = apply(X, 2, function(x){ sum(x^2, na.rm = TRUE) })
      s0[1] = (1 / sum(xx) ) * R2 * var(y, na.rm = TRUE) * ( (nu0[1] + 2) / nu0[1])
      s0[2] = (1 / sum(xx) ) * R2 * var(y, na.rm = TRUE) * ( (nu0[2] + 2) / nu0[2])

   } else {
     stop("Error")
   }  

   return(s0)
}   


   
