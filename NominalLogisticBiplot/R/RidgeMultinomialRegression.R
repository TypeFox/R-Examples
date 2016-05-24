# file NominalLogisticBiplot/R/RidgeMultinomialRegression.R
# copyright (C) 2012-2013 J.L. Vicente-Villardon and J.C. Hernandez
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#Function that calculates an object with the fitted multinomial logistic regression for a nominal variable.
#It compares with the null model, so that we will be able to compare which model fits better the variable.
#----------------------Parameters--------------
  #y: response nominal variable.
  #x: matrix with independent variables.
  #penalization: value to correct the separation problem through the ridge regression
  #cte : it will be true if the model has an independent term.
  #tol : value to decide if the algorith should continue
  #maxiter : value to decide if the algorith should continue
  #showIter: boolean parameter if we want to see values in each iteration of the process.
RidgeMultinomialRegression <- function(y, x, penalization = 0.2, cte = TRUE, tol = 1e-04, maxiter = 200,showIter=FALSE){
  if(is.matrix(x)){
  	n <- nrow(x)
  }else{
    n <- length(x)
  }

	Model=polylogist(y, x, penalization=penalization, tol=tol, maxiter=maxiter, show=showIter)
	Null=polylogist(y, matrix(1,n,1), penalization=penalization, tol=tol, maxiter=maxiter, show=showIter, cte = FALSE)

	Fit=Model
	Fit$NullDeviance=Null$Deviance
	Dif=(Null$Deviance - Model$Deviance)
	Fit$Difference=Dif
	Fit$df=length(Model$beta)-length(Null$beta)
	Fit$p=1-pchisq(Dif, df =  Fit$df)
	Fit$CoxSnell=1-exp(-1*Dif/n)
  Fit$Nagelkerke=Fit$CoxSnell/(1-exp((Null$Deviance/(-2)))^(2/n))
	Fit$MacFaden=1-(Model$Deviance/Null$Deviance)
	     predicted=array(0,n)
	     for(i in 1:n){       
	       predicted[i] = sum(which(Model$fitted[i,]==max(Model$fitted[i,]))==y[i])
       }
  Fit$PercentCorrect = sum(predicted)/n *100
	return(Fit)
}

