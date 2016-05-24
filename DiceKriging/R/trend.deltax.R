trend.deltax <- function(x, model, h=sqrt(.Machine$double.eps)){
  
  n <- model@n
  d <- model@d
  
  x <- matrix(x, nrow=1)
  if (length(x)!=d) {
    stop("x must be a vector")
  }
  names.x <- colnames(model@X)
  colnames(x) <- names.x
  rownames(x) <- NULL

  formula <- model@trend.formula
  
  if (formula==~1){ # OK case
    grad.intercept <- matrix(0, nrow=1, ncol=d,
                             dimnames=list("(Intercept)", 1:d))
    return(grad.intercept)      
  } 
  
  formula.linear <- drop.response(~., data=data.frame(x))
  formula.quad <- drop.response(~.^2, data=data.frame(x))
  
  if ((formula==formula.linear) | (formula==formula.quad)) {
    grad.intercept <- matrix(0, nrow=1, ncol=d,
                             dimnames=list("(Intercept)", 1:d))
    grad.linear <- diag(d)
    rownames(grad.linear) <- names.x
    grad.linear <- rbind(grad.intercept, grad.linear)
    if (formula==formula.linear) {
      return(grad.linear)
    }  
    grad.inter <- matrix(0, nrow=as.integer(d*(d-1)/2), ncol=d)
    names.f <- colnames(model.matrix(~.^2, data=data.frame(x)))
    names.inter <- names.f[-(1:(d+1))]
    for (j in 1:d){
      index <- grep(names.x[j], names.inter)
      grad.inter[index,j] <- x[-j]
    }
    rownames(grad.inter) <- names.inter
    return(rbind(grad.linear, grad.inter))
  } # end analytic cases
  
  A <- matrix(x, nrow=d, ncol=d, byrow=TRUE)
  colnames(A) <- colnames(x)
  Apos <- A+h*diag(d)
  Aneg <- A-h*diag(d)
  newpoints <- data.frame(rbind(Apos, Aneg))
  f.newdata <- model.matrix(model@trend.formula, data = newpoints)
  f.deltax <- (f.newdata[1:d,]-f.newdata[(d+1):(2*d),])/(2*h)
  f.deltax <- t(f.deltax) 
  
  return(f.deltax)

}