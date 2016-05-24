# file NominalLogisticBiplot/R/NominalLogBiplotEM.R
# copyright (C) 2012-2013 J.C. Hernandez and J.L. Vicente-Villardon
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

#This function is an alternated algorith for calculating en two steps the row coordinates,
# or also known as ability (E-step), and parameters of each model for the variables (M-step).
# The algorith uses multidimensional Gauss-Hermite quadrature and Ridge Regression to solve the
# separation problem in logistic regression.
NominalLogBiplotEM <- function(x, dim = 2, nnodos = 10, tol = 1e-04, maxiter = 100, penalization = 0.2,
                               initial=1,alfa=1,Plot=FALSE,showResults=FALSE) {
  initials = c("MCA","MIRT")
  if(is.numeric(initial)){
    initial = initials[initial]                                                                                         
  }

	Q = multiquad(nnodos, dim)
	X = Q$X
	A = Q$A
	q = dim(X)[1]
  p = dim(x)[2]
  Ncats=ColMax(x)
  Maxcat=max(Ncats)
	G = NominalMatrix2Binary(x)
  s = dim(G)[1]
	n = dim(G)[2]
	VariableModels = 0

  par = array(0, c(Maxcat-1, dim + 1, p))    

	if(showResults){print("Calculating initial coordinates")}
  if(initial == "MCA"){
      dimension = dim
     	corr = afc(G, dim = dimension,alpha=alfa)
    	ability = corr$RowCoordinates[, 1:dim]
  }
  if(initial == "MIRT"){
      technical = list()
      technical$NCYCLES = maxiter            
      nominalHighLow = matrix(0,2,p)
      for(i in 1:p){
          nominalHighLow[1,i] = as.numeric(max(x[,i]))
          nominalHighLow[2,i] = as.numeric(min(x[,i]))
      }                                                                           
      mod = mirt(x ,dim,"nominal",nominal.highlow=nominalHighLow,technical=technical)
      ability = fscores(mod,full.scores=T)[,1:dim] 
  }  
	if(showResults){print("Initial coordinates calculated")}
	
	logLikold=0
	for (j in 1:p){
	  if(showResults){print(paste("Adjusting variable ",j,sep=""))} 
  	model = RidgeMultinomialRegression(x[, j], ability, penalization = penalization, tol = tol, maxiter = maxiter,showIter = showResults)  	
  	vmodel = model
  	vmodel$name = dimnames(x)[[2]][j]
   	VariableModels=cbind(VariableModels,vmodel)
  	par[1:(Ncats[j]-1),,j] =vmodel$beta
  	logLikold=logLikold+vmodel$logLik
  }
  VariableModels=VariableModels[,2:(p+1)]
  
	if (Plot){
  	plot(ability[,1],ability[,2])
    text(ability[,1],ability[,2],1:q, pos=4, offset=0.2)
  }
  
	error = 1
	iter = 0      
  if(showResults) print("Starting iterations in the algorithm") 
	while ((error > tol) & (iter < maxiter)) {
		# E-step
		iter = iter + 1
		PT= EvalPolylogist(X, par, Ncats)
		L = matrix(1, s, q)
		for (l in 1:s)
     for (k in 1:q)
       L[l, k] = prod(PT[k,]^G[l,])
		Pl = L %*% A
		ability = matrix(0, s, dim)
		for (l in 1:s)
     for (j in 1:dim){
			for (k in 1:q)
         ability[l, j] = ability[l, j] + X[k, j] * L[l, k] * A[k]
			ability[l, j] = ability[l, j]/Pl[l]
		}
  	###################### M-step  
  	logLik=0
    for (j in 1:p){
  	  if(showResults){print(paste("Adjusting variable ",j,sep=""))}     
 	  	model = RidgeMultinomialRegression(x[, j], ability, penalization = penalization, tol = tol, maxiter = maxiter,showIter = showResults)  	
 	   	vmodel = model
    	vmodel$name = dimnames(x)[[2]][j]  
 	  	VariableModels[,j] = vmodel
    	par[1:(Ncats[j]-1),,j] =vmodel$beta
    	logLik=logLik+vmodel$logLik
    }
		error = abs((logLik - logLikold)/logLik)
		logLikold = logLik
		if(showResults){
       print(paste("Iteration ",iter,"- Log-Lik:",logLik," - Change:",error,sep=""))		
    }
	}
  print(paste("Iteration ",iter,"- Log-Lik:",logLik," - Change:",error,sep=""))
	if (Plot)
    	for (j in 1:p){
      	dev.new()
      	plot(ability[,1],ability[,2],col=x[,j])
        text(ability[,1],ability[,2],1:q,col=x[,j])
      }
	d = sqrt(rowSums(par^2) + 1)
	model = list()
	model$RowCoordinates = ability
	model$ColumnModels = VariableModels
	class(model) = "nominal.logistic.biplot.EM"
	return(model)
}

