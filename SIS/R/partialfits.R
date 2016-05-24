obtain.beta<-function(x, y, ix0, family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic){  
if(tune == "cv"){
   if(penalty == "lasso"){
      cv.fit = cv.glmnet(x[,ix0], y, family=family, type.measure=type.measure, nfolds=nfolds)
      if(family == "cox")
         coef.beta = coef(cv.fit,s="lambda.min")                         # extract coefficients at a single value of lambda, there is no intercept in the model
      else
         coef.beta = coef(cv.fit,s="lambda.min")[-1]                     # extract coefficients at a single value of lambda, delets the intercept   
   }else{
       X = as.matrix(x[,ix0])       
       cv.fit = cv.ncvreg(X=X, y=y, family=family, penalty=penalty, gamma=concavity.parameter, nfolds=nfolds)
       coef.beta = cv.fit$fit$beta[-1,cv.fit$min]   
   }                        
}else{
    n = nrow(x)
    if(penalty == "lasso"){
       reg.fit = glmnet(x[,ix0], y, family=family)
       coef.beta =  reg.fit$beta
       dev = deviance(reg.fit)
       reg.df = reg.fit$df
     }
     else{
        reg.fit = ncvreg(as.matrix(x[,ix0]), y, family=family, penalty=penalty, gamma=concavity.parameter)
        coef.beta = as.matrix(reg.fit$beta[-1,])
        dev = loglik(as.matrix(x[,ix0]), y, reg.fit$beta, family=family)        
        reg.df = getdf(coef.beta)
     }
     if(tune == "aic") {obj = dev + 2*reg.df}
     if(tune == "bic") {obj = dev + log(n)*reg.df}
     if(tune == "ebic"){obj = dev + log(n)*reg.df + 2*gamma.ebic*log(choose(length(ix0), reg.df))}
     ind = which.min(obj)
     if (ncol(coef.beta) > 1) coef.beta =  coef.beta[,ind]
     else coef.beta =  coef.beta[ind]
}  
return(as.vector(coef.beta))
}   

final.fit<-function(x, y, family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic){  
if(tune == "cv"){
   if(penalty == "lasso"){
      fit = cv.glmnet(x, y, family=family, type.measure=type.measure, nfolds=nfolds)
      ind = which(fit$lambda==fit$lambda.min)
   }else{
       X = as.matrix(x)       
       fit = cv.ncvreg(X=X, y=y, family=family, penalty=penalty, gamma=concavity.parameter, nfolds=nfolds)  
       ind = fit$min
   }   
                         
}else{
    n = nrow(x)
    if(penalty == "lasso"){
       fit = glmnet(x, y, family=family)
       dev = deviance(fit)
       reg.df = fit$df
     }
     else{
        fit = ncvreg(as.matrix(x), y, family=family, penalty=penalty, gamma=concavity.parameter)
        coef.beta = as.matrix(fit$beta[-1,])
        dev = loglik(as.matrix(x), y, fit$beta, family=family)        
        reg.df = getdf(coef.beta)
     }
     if(tune == "aic") {obj = dev + 2*reg.df}
     if(tune == "bic") {obj = dev + log(n)*reg.df}
     if(tune == "ebic"){obj = dev + log(n)*reg.df + 2*gamma.ebic*log(choose(dim(as.matrix(x))[2], reg.df))}
     ind = which.min(obj)
}  
return(list(fit=fit, ind=ind))
}  
  
loglik=function(X, y, beta, family){
  K = dim(beta)[2]
  link = cbind(1,X)%*%beta
  yrep = repmat(y, 1, K) 
  if(family=="gaussian") return(apply((yrep-link)^2, 2, sum))
  if(family=="poisson") return(apply(exp(link)-yrep*link, 2, sum))
  if(family=="binomial") return(apply(log(1+exp(link))-yrep*link, 2, sum))
}

repmat = function(X,m,n){
##R equivalent of repmat (matlab)
X = as.matrix(X)
mx = dim(X)[1]
nx = dim(X)[2]
matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

getdf = function(coef.beta){
  apply(abs(coef.beta)>1e-10,2, sum)
} 
