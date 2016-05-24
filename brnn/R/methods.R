# file brnn/methods.R
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

#A package for Bayesian Regularized Neural Networks
#Author: Paulino Perez Rodriguez
#Madison, WI, Sep. 2012
#Birmingaham, Alabama, Jan. 2013

print.brnn=function(x,...)
{
  	if(!inherits(x, "brnn")) stop("This function only works for objects of class `brnn'\n");
  	cat("A Bayesian regularized neural network \n");
  	cat(paste(x$p,"-",x$neurons,"- 1 with",x$npar, "weights, biases and connection strengths\n",sep=" "));
  	cat("Inputs and output were", ifelse(x$normalize,"","NOT"),"normalized\n",sep=" ");
  	cat("Training finished because ",x$reason,"\n");
}

summary.brnn=function(object,...)
{
	object
}

predict.brnn=function(object,newdata,...)
{
   y=NULL;

   if(!inherits(object,"brnn")) stop("This function only works for objects of class `brnn \n'");
   
   if (missing(newdata) || is.null(newdata)) 
   {
        
        y=predictions.nn.C(vecX=as.vector(object$x_normalized),n=object$n,p=object$p,
                           theta=object$theta,neurons=object$neurons,cores=1);
        
        #predictions in the original scale
        if(object$normalize)
        {
            y=un_normalize(y,object$y_base,object$y_spread)
        }
        
   }else{

	if(inherits(object,"brnn.formula"))
	{
		 ## The model was fitted using formula interface
		 newdata = as.data.frame(newdata)

            	 ## Obtain the design matrix for the prediction problem
            	 Terms = delete.response(object$terms)
            	 mf=model.frame(Terms, newdata, na.action = na.omit, xlev = object$xlevels)
            	 if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, mf)
            	 x = model.matrix(Terms, mf, contrasts = object$contrasts)
            	 xint = match("(Intercept)", colnames(x), nomatch=0L)
            	 if(xint > 0L) x = x[, -xint, drop=FALSE] 
	}else{
		# Matrix fit
		if(is.null(dim(newdata))) dim(newdata) = c(1L, length(newdata)) # a row vector
        	x = as.matrix(newdata)     # to cope with dataframes
        	if(any(is.na(x))) stop("missing values not allowed in 'newdata' \n")
	}
        
	if(ncol(x)!=object$p) stop("The number of predictors used to fit the model and those in `newdata' does not match\n");
        
        #normalize using the information in the trainig set
        if(object$normalize)
        {
        	x=normalize(x,base=object$x_base,spread=object$x_spread)
        }
        
        y=predictions.nn.C(vecX=as.vector(x),n=nrow(x),p=ncol(x),
                         theta=object$theta,neurons=object$neurons,cores=1);                 
        if(object$normalize)
        {
        	y=un_normalize(y,object$y_base,object$y_spread)
        }
   }
   return(y)
}

#FIXME:
#Return the output
coef.brnn=function(object, ...)
{
	if(!inherits(object, "brnn")) stop("This function only works for objects of class `brnn'\n");
        theta=object$theta
        theta_names=as.character()        

	#Loop for getting weights, biases and connection strengths for each neuron
	for(i in 1:(object$neurons))
    	{
           theta_names=c(theta_names,paste("w[",i,"]",sep=""),paste("b[",i,"]",sep=""),paste(paste("beta[",i,",",sep=""),1:object$p,"]",sep=""))
    	}
        theta=as.vector(unlist(theta))
        names(theta)=theta_names
        theta
}

predict.brnn_extended=function(object, newdata,...)
{
	if(!inherits(object, "brnn_extended")) stop("This function only works for objects of class `brnn.extended'\n");
	y=NULL;

	if (missing(newdata) || is.null(newdata))
   	{

        	y=predictions.nn.C(vecX=as.vector(object$x_normalized),n=object$n,p=object$p,theta=object$theta1,neurons=object$neurons1,cores=1) +
                  predictions.nn.C(vecX=as.vector(object$z_normalized),n=object$n,p=object$q,theta=object$theta2,neurons=object$neurons2,cores=1)

        	#predictions in the original scale
        	if(object$normalize)
        	{
            		y=un_normalize(y,object$y_base,object$y_spread)
        	}

   }else{
	
	if(inherits(object,"brnn_extended.formula"))
        {
            Terms = delete.response(object$mtx)
            mf=model.frame(Terms, newdata, na.action = na.omit)
            x = model.matrix(Terms, mf,contrasts = object$contrastsx,xlev = object$xlevels)
            if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, mf)
            xint = match("(Intercept)", colnames(x), nomatch=0L)
            if(xint > 0L) x = x[, -xint, drop=FALSE]
            
            Terms = delete.response(object$mtz)
            mf=model.frame(Terms, newdata, na.action = na.omit)
            z = model.matrix(Terms, mf,contrasts = object$contrastsz,xlev = object$zlevels)
            if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, mf)
            zint = match("(Intercept)", colnames(z), nomatch=0L)
            if(zint > 0L) z = z[, -zint, drop=FALSE]
            
            xz=as.matrix(cbind(x,z))
            
        }else{
		#Matrix fit
        	if(is.null(dim(newdata))) dim(newdata) = c(1L, length(newdata)) # a row vector
        	xz = as.matrix(newdata)     # to cope with dataframes
        	if(any(is.na(xz))) stop("missing values not allowed in 'newdata' \n")
	}

        if(ncol(xz)!=(object$p+object$q)) stop("The number of predictors used to fit the model and those in `newdata' does not match\n");
        
        #normalize using the information in the trainig set
        if(object$normalize)
        {
                x=normalize(xz[,c(1:object$p)],base=object$x_base,spread=object$x_spread)
                z=normalize(xz[,c((object$p+1):(object$p+object$q))],base=object$z_base,spread=object$x_spread)
        }

         y=predictions.nn.C(vecX=as.vector(x),n=nrow(xz),p=object$p,theta=object$theta1,neurons=object$neurons1,cores=1) +
           predictions.nn.C(vecX=as.vector(z),n=nrow(xz),p=object$q,theta=object$theta2,neurons=object$neurons2,cores=1)

        if(object$normalize)
        {
                y=un_normalize(y,object$y_base,object$y_spread)
        }
   }
   return(y)	
}

print.brnn_extended=function(x,...)
{
	if (!inherits(x, "brnn_extended")) stop("This function only works for objects of class `brnn_extended \n'");
  	cat("A Bayesian regularized neural network \n");
  	cat(paste(x$p,"-",x$neurons1,"-",x$neurons2,"- 1 with",x$npar1+x$npar2, "weights, biases and connection strengths\n",sep=" "));
  	cat("Inputs and output were", ifelse(x$normalize,"","NOT"),"normalized\n",sep=" ");
  	cat("Training finished because ",x$reason,"\n");
}

summary.brnn_extended=function(object,...)
{
	object
}

coef.brnn_extended=function(object,...)
{
	if(!inherits(object, "brnn_extended")) stop("This function only works for objects of class `brnn_extended'\n");
        theta1=object$theta1
        theta1_names=as.character()

        #Loop for getting weights, biases and connection strengths for each neuron
        for(i in 1:(object$neurons1))
        {
           theta1_names=c(theta1_names,paste("w[",i,"]",sep=""),paste("b[",i,"]",sep=""),paste(paste("beta[",i,",",sep=""),1:object$p,"].x",sep=""))
        }
        theta1=as.vector(unlist(theta1))
        names(theta1)=theta1_names

        theta2=object$theta2
        theta2_names=as.character()
        for(i in 1:(object$neurons2))
        {
           theta2_names=c(theta2_names,paste("w[",i,"]",sep=""),paste("b[",i,"]",sep=""),paste(paste("beta[",i,",",sep=""),1:object$q,"].z",sep=""))
        }
        theta2=as.vector(unlist(theta2))
        names(theta2)=theta2_names
        out=list()
        out[[1]]=theta1
        out[[2]]=theta2
        out
}
