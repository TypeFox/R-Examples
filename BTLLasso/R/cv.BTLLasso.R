
cv.BTLLasso <- function(Y, X, folds = 10, lambda, control = BTLLasso.ctrl(), cores = folds,
                        trace = TRUE, trace.cv = TRUE){

  
  vardiffs <- abs(apply(X,2,var)-1)
  
  if(!is.matrix(X))
    stop("X has to be a matrix")
  
  if(!is.matrix(Y))
    stop("Y has to be a matrix")
  
  ### extract all control arguments
  adaptive <- control$adaptive
  norm <- control$norm
  epsilon <- control$epsilon
  lambda2 <- control$lambda2
  c <- control$c
  penal.diffs <- control$penal.diffs
  return.design <- control$return.design
  ###
  
  
  if(trace.cv){cat("Full model","\n")}
  m.all <- BTLLasso(Y=Y, X=X, lambda=lambda, control=BTLLasso.ctrl(adaptive = adaptive, 
              norm = norm, epsilon = epsilon, lambda2 = lambda2, c = c, penal.diffs = penal.diffs, 
              return.design = TRUE), trace = trace)
 
  design <- m.all$design
  Y.ord <- m.all$Y
  q <- m.all$q
  k <- q+1
  
  n <- nrow(X)
  I <- ncol(Y)
  m <- (1 + sqrt(1+8*I))/2
  subjects <- 1:n
  
  
  n.cv <- rep(floor(n/folds),folds)
  rest <- n%%folds
  if(rest>0){
    n.cv[1:rest] <- n.cv[1:rest] +1}
  
  which.fold <- rep(1:folds,n.cv)
  
  id.fold <- sample(which.fold,n,replace=FALSE)
  
  
  
  
  cv.fun <- function(ff){
  
    if(trace.cv){cat("CV-fold:",ff,"out of",folds,"\n")}
    
    X.train <- X[-which(id.fold==ff),,drop=FALSE]
    
    Y.train <- Y[-which(id.fold==ff),,drop=FALSE]
    
    X.test <- X[which(id.fold==ff),,drop=FALSE]
    
    Y.test <- c(Y[which(id.fold==ff),,drop=FALSE])
    
    m.fold <- BTLLasso(Y=Y.train, X=X.train, lambda=lambda, control = BTLLasso.ctrl(adaptive = adaptive, 
                          norm = norm, epsilon = epsilon, lambda2 = lambda2, c = c,
                          penal.diffs = penal.diffs, return.design = FALSE), trace = trace)

    coef.fold <- m.fold$coefs

    preds <- c()
    for(u in 1:length(lambda)){
      preds <- cbind(preds,predict.BTLLasso(X.test,coef.fold[u,],q,m))
    }
    
    yhelp <- rep(Y.test,each=k)
    yhelp <- as.numeric(yhelp==rep(1:k,length(Y.test)))
    
    deviance <- - 2*colSums(yhelp*log(preds))
    deviance
  }
    
  if(cores>1){
  cl <- makeCluster(cores,outfile="")
  
  clusterExport(cl, varlist = c("X","Y","id.fold","lambda","adaptive","norm","epsilon",
                                "lambda2","c","penal.diffs") 
                ,envir = sys.frame(sys.nframe()))
                  
  deviances <- rowSums(parSapply(cl, seq(folds), cv.fun))
  stopCluster(cl)
  }else{
    deviances <- rowSums(sapply(seq(folds), cv.fun))
  }
  

  
  ret.list <- list(coefs = m.all$coefs, coefs.repar = m.all$coefs.repar, logLik =m.all$logLik, 
                   design = m.all$design, Y = Y, q = q, acoefs = m.all$acoefs, 
                   response = m.all$response, n = n, I = I, m = m, p = m.all$p, X = X,  
                   n.theta = m.all$n.theta, lambda = lambda, deviances = deviances, folds = folds,
                    labels = m.all$labels, epsilon = epsilon)
  
  class(ret.list) <- c("cv.BTLLasso","BTLLasso")
  
  return(ret.list)
}
