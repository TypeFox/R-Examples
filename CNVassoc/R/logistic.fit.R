logistic.fit<-function(x,y,weights,tol=1e-8,max.iter=25,verbose=FALSE){

  if (!is.matrix(x))
    x <- as.matrix(x)

  nvars <- NCOL(x)
  nobs <- NROW(x)

  mu <- drop((weights * y + 0.5)/(weights + 1))
  eta <- log(mu/(1-mu))

  iter<-0  
  prev.dev <- dev.increase <- Inf
  
  while(dev.increase>=tol & iter<=max.iter){  
    
    iter <- iter+1

    z <- eta+(y-mu)/(mu*(1-mu))

    ww <- sqrt(weights*mu*(1-mu))
    
    fit <- .Fortran("dqrls", qr = x * ww, n = nobs, 
            p = nvars, y = ww * z, ny = as.integer(1), tol = min(1e-07, 
            glm.control()$epsilon/1000), coefficients = double(nvars), 
            residuals = double(nobs), effects = double(nobs), 
            rank = integer(1), pivot = 1:nvars, qraux = double(nvars), 
            work = double(2 * nvars), PACKAGE = "CNVassoc")
    
    coeffic <- fit$coefficients
    eta <- drop(x%*%coeffic)
    mu <- 1/(1+exp(-eta))
    
    cur.dev <- -2*sum(y*eta+log(1-mu))
    dev.increase <- abs(prev.dev-cur.dev)
    prev.dev <- cur.dev
    
    if (verbose) 
      cat("Deviance iter",iter,"=",cur.dev,"\n")

  } 

  coeffic<-as.vector(coeffic)
  names(coeffic)<-colnames(x)
  
  return(coeffic)

}
