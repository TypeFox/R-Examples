predict.RAMP <- function(object, newdata=NULL, type=c("link", "response", "class"), allpath=FALSE,...){
# object: a fitted object of class inheriting from "RAMP".
# X: optionally, a data matrix with which to predict. If omitted, the data in object will be used.
# type: the type of prediction required. The default is on the scale of the linear predictors; the #alternative "response" is on the scale of the response variable. Thus for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives the predicted probabilities.  type="class" gives the class lable for binomail distribution.
    X = newdata
    oldX = object$X
    if(is.null(X)) X=oldX   
    n = nrow(X)
    p = ncol(X)
	mainind.list = object$mainInd.list
	interind.list = object$interInd.list
	a0.list = object$a0
	k = length(mainind.list)-1 ###total number of models visited
	
	eta = matrix(0,n,k)	
	model.list = as.list(NULL)
	count = 0
	for(i in 1:k){
		#cat('i= ', i, ' out of ', k,'\n')
		mainind = mainind.list[[i]]
		Xmain = X[,mainind]
		count = count+1
		Xinter = NULL
		Xi = Xmain
	    coef = object$beta.m.mat[mainind,i]
		if(i<=length(interind.list)){
		interind = interind.list[[i]]		
        if(length(interind)>0&&(!(length(interind)==1&&interind=='None'))){##candidate interaction terms
             	    for(indInter in 1:length(interind)){
             	    pair = as.numeric(strsplit(interind[indInter],'X')[[1]][2:3])
             	    Xinter = cbind(Xinter, X[,pair[1]]*X[,pair[2]])
	              	}
	              	
	    
	    if(!is.na(pair[1])&&!is.na(pair[2])){
	    Xi = cbind(Xmain,Xinter)   
	    coef = c(object$beta.m.mat[mainind,i],object$beta.i.mat[[i]])
		} 
		}
		}
		eta[,i] = a0.list[i]
		if(length(mainind)>0){
	    eta[,i] = a0.list[i]+as.matrix(Xi)%*%coef
	    }
		}
	if(allpath==FALSE){
	eta = eta[,object$cri.loc]
	}
	if(match.arg(type)=="response") return(exp(eta)/(1+exp(eta)))
    if(match.arg(type)=="link") return(eta)
    if(match.arg(type)=="class") return(as.integer(eta>0))

}
