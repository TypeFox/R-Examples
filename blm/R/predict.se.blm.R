predict.se.blm <- function(object,newdata,se=FALSE){
  
  if(missing(newdata))
    newdata <- object@data
    
  if(class(object)=="lexpit"){
    X <- model.matrix(update(object@formula.linear,NULL~.),newdata)
    Z <- model.matrix(update(object@formula.expit,NULL~.),newdata)
    
    x.has.intercept <- attr(terms(object@formula.linear),"intercept")==1
    
    if(x.has.intercept) X <- X[,-1]
    
    if(!is.matrix(X)) X <- cbind(X)
    if(!is.matrix(Z)) Z <- cbind(Z)
    
    p <- X%*%object@coef.linear+expit(Z%*%object@coef.expit)
  }
  else{
    X <- model.matrix(update(object@formula,NULL~.),newdata)
    p <- X%*%object@coef
  }
  
  if(!se){
    return(as.numeric(p))
  }
  else{
    
    if(class(object)=="blm"){
      se <- sqrt(apply(X,1,function(x) t(x)%*%object@vcov%*%x))
    }
    else{
      var <- apply(X,1,function(x) t(x)%*%object@vcov.linear%*%x)      
      d.expit <- function(x) expit(x)*(1-expit(x))
      g <- function(x) x*d.expit(x%*%object@coef.expit)
      var.expit <- apply(Z,1,function(x) t(g(x))%*%object@vcov.expit%*%g(x))
      se <- sqrt(var+var.expit)
    }
    
    return(data.frame(fit = p, se = se))
  }
}

setMethod("predict","blm",predict.se.blm)
setMethod("predict","lexpit",predict.se.blm)
