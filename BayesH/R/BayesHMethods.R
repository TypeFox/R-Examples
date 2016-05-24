######################################################################
#
# BayesHMethods.R 
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
# Contains: predict.BayesH, plot.BayesH and summary.BayesH
# Remarks: This code is an adaptation of code written by Perez and    
# De Los Campos (2014).       
#
#######################################################################
predict.BayesH =function(object,newdata, ...)
{
    if (!inherits(object, "BayesH")){
       stop("This function only works for objects of class `BayesH'")
    }
    
    if(missing(newdata) || is.null(newdata)) {
	return(object$fitted.values)
    } else if(!is.null(newdata) && (is.matrix(newdata) || is.data.frame(newdata))){
        return(as.matrix(newdata)%*%object$coef)
    } else {
       stop("You must define the newdata as matrix or data.frame")
    }    
}


plot.BayesH =function(x,...)
{
  ### Check that object is compatible
  if(!inherits(x, "BayesH")){
     stop("This function only works for objects of class `BayesH'")
  }
  #limits<-range(c(x$y,x$fitted.values),na.rm=TRUE)
  plot(x$y,x$fitted.values,main="Training",xlab='Response',ylab='Prediction', pch=16)
  abline(a=0,b=1,lty=3)
}


summary.BayesH = function(object,...)
{
    if(!inherits(object, "BayesH")) {
        stop("This function only works for objects of class `BayesH'")
    }
    tmp = paste('-------------------Summary of model -------------------------')
    cat(tmp,'\n\n')
    cat('Variance of observed response (TRN)=', round(var(object$y,na.rm=TRUE),4),'\n') 
    cat('Error variance=',round(mean(object$sigmaGibbs, na.rm=TRUE),4),'\n')
    n = length(object$y)
    if(any(is.na(object$y)))
    {     
     		tst = which(is.na(object$y))
     		cat('N-TRN=',n-length(tst), 'N-TST=',length(tst),'\n')     
     		cat('Correlation TRN=',round(cor(object$y,object$fitted.values, use = "complete.obs"),4),'\n')
    } else{
       cat('N-TRN=',n,'N-TST=0', '\n\n')
    }  
}
