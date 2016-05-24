tune.fit<-function(x, y, family = c("gaussian","binomial","poisson","cox"), penalty=c("SCAD","MCP","lasso"), concavity.parameter = switch(penalty, SCAD=3.7, 3), tune = c("cv","aic","bic","ebic"), nfolds = 10, 
                   type.measure = c("deviance","class","auc","mse","mae"), gamma.ebic = 1){  
                           
if(is.null(x)||is.null(y)) stop("The data is missing!")
        
this.call=match.call()
family = match.arg(family)
penalty = match.arg(penalty)  
if(class(concavity.parameter) != "numeric") stop("concavity.parameter must be numeric!")
tune = match.arg(tune)   
if(class(nfolds) != "numeric") stop("nfolds must be numeric!")
type.measure = match.arg(type.measure)                          
                           
if(tune == "cv"){
   if(penalty == "lasso"){
      cv.fit = cv.glmnet(x, y, family=family, type.measure=type.measure, nfolds=nfolds)
      if(family == "cox")
         coef.beta = coef(cv.fit,s="lambda.min")                         
      else
         coef.beta = coef(cv.fit,s="lambda.min")                       # extract coefficients at a single value of lambda, including the intercept
   }else{
       cv.fit = cv.ncvreg(x, y, family=family, penalty=penalty, gamma=concavity.parameter, nfolds=nfolds)
       coef.beta = cv.fit$fit$beta[,cv.fit$min]                        # extract coefficients at a single value of lambda, including the intercept
   }                        
}else{
    n =  nrow(x)
    if(penalty == "lasso"){
       reg.fit = glmnet(x, y, family=family)
       coef.beta = reg.fit$beta                                        # extract coefficients at all values of lambda, NOT including the intercept
       a0 = reg.fit$a0                                                 # extract intercept at all values of lambda                      
       dev = deviance(reg.fit)
       reg.df = reg.fit$df
     }
     else{
        reg.fit = ncvreg(x, y, family=family, penalty=penalty, gamma=concavity.parameter)
        coef.beta =  reg.fit$beta                                      # extract coefficients at all values of lambda, including the intercept
        dev = loglik(x, y, reg.fit$beta, family=family)
        reg.df = getdf(coef.beta)
     }
     if(tune == "aic") {obj = dev + 2*reg.df}
     if(tune == "bic") {obj = dev + log(n)*reg.df}
     if(tune == "ebic"){obj = dev + log(n)*reg.df + 2*gamma.ebic*log(choose(dim(x)[2], reg.df))}
     ind = which.min(obj)
     coef.beta = coef.beta[,ind]
     if(penalty == "lasso"){
        a0 = a0[ind] 
        coef.beta = c(a0, coef.beta)                                   # combines intercept and coeffcients at all values of lambda
     }      
} 
ix = sort(which(coef.beta[-1]!=0))
fit = coef.beta[-1][ix]
a0 = coef.beta[1]
return(list(ix = ix, a0=a0, fit=fit))
} 
