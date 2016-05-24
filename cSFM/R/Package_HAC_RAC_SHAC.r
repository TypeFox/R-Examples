# functions used to estimate the model 
# given an input data set, output the estimation results


# library(sn)
# source("BasicSetting.r")

# these are particularly for case=2, where both the mean and variance are covariant dependent; 
#                            case=1: all functional parameters are covariant independent.



### kpbb: generate smooth matrices, output the smooth matrices for all the parameter functions : mu, logvar, gamma(trasformed skewness)
### use regression spline with degree = degree.poly.                  
# input:    timepoints - tp(must be sorted);
#           covariate - cp (must be sorted); 
#           case number - specifying the dependence on covariate; 1 for covariant independent, and 2 for covariant dependent in mean and variance
#           degree.poly - degree of polynomial, default as 3 (cubic spline); same for different parameters.
# output: list(sm.mu, sm.logvar, sm.gamma, index) - 3 smooth matrics and a index vector with 3 distint values: 1 for mu, 2 for logvar, 3 for gamma (skewness) 


# Kronecker product B-spline Basis

kpbb <- function(tp, cp = NULL, nknots.tp, 
                 nknots.cp = NULL, sub.case = 2, degree.poly = 3){
  # require(splines)
  if(is.unsorted(tp)) {stop("timepoints are not sorted")}     
  basis.time <- bs(tp, df= degree.poly + nknots.tp + 1, degree= degree.poly, intercept=T) 
  
  if (sub.case == 1)
  {
    basis <- basis.time
  }
  
  if (sub.case == 2){  #cov dependent
    if (is.null(cp)){stop("covariant not provided")}
    if (is.unsorted(cp)) { stop("covariant not sorted")}    
    if (is.null(nknots.cp)) { 
      stop("number of knots for covariate is not specified")
    }       
    
    basis.cov  <- bs(cp, df= degree.poly + nknots.cp + 1, degree= degree.poly, intercept=T)     
    basis.tensor <- kronecker(basis.time, basis.cov)    
    
    new.attr <- attributes(basis.tensor)
    new.attr$tp <- attributes(basis.time)
    new.attr$cp <- attributes(basis.cov)
    attributes(basis.tensor) <- new.attr
    
    basis <- basis.tensor
  }   
  
  class(basis) <- c(class(basis), "kpbb")
  return(basis)
}

  



### case2.unmll.optim: objective function to be minimize, which is the -loglikelihood in terms of beta's. 
#INPUT: beta - smoothing coefficients as a vector; 
#      dataone - observation as a matrix n by m (n: number of subjects; m: number of timepoints);
#      Basis.list - a list of 3 components; each one corresponds to a smooth matrix (for mu, logvar, and gamma); 
#      cate - category of model to be considered; 1 for full model, 2 for the model when the shape is fixed at 0; 
#OUTPUT: -loglikehood at beta when data = dataone and the Basis.list is used. 
        
  case2.unmll.optim <- function(beta, dataone, Basis.list, cate=1){ 
    y <- dataone
    n <- nrow(y)
    m <- ncol(y)    
    # index: indicate the length and type of beta's 
    #   1 for beta_mu, 2 for beta_logvar, 3 for beta_gamma
    index = unlist(lapply(c(1:3), function(k) rep(k, ncol(Basis.list[[k]]))))  
    beta.list = lapply(c(1:3), function(k) beta[index == k])
    
    if (cate == 2){
      beta.list[[3]] <- rep(0,length(beta[index == 3]))}
    
    if (cate == 3){
      beta.list[[1]] <- rep(0,length(beta[index == 1]))}
    
    est <- array(NA, c(3,n,m))              
    #estimation of B %*% beta; NOT var and skewness but logvar and gamma 
    est[1,,]  <- matrix(Basis.list[[1]] %*% beta.list[[1]], n, m, byrow=F)
    est[2,,] <-  matrix(Basis.list[[2]] %*% beta.list[[2]], n, m, byrow=F)
    est[3,,] <-  matrix(c(Basis.list[[3]] %*% beta.list[[3]]), n, m, byrow=T)               
    
    mu <- est[1,,]
    sigma <- exp(est[2,,]/2)
    x <- (dataone - mu)/sigma
    gamma1 <- (pnorm(est[3,,])*2 - 1)*0.99527174643
    
    unmll <- -sum(g((y- mu)/sigma,gamma1, log = TRUE) - est[2,,]/2)    #unmll
    
    return(unmll)
  }

  # gradient of case2.unmll, same input as case.unmll.optim but return the gradient
  
  case2.gr <- function(beta, dataone, Basis.list, cate = 1){
    y <- dataone
    n <- nrow(y)
    m <- ncol(y)
    index = unlist(lapply(c(1:3), function(k) rep(k, ncol(Basis.list[[k]]))))  
    beta.mu <- beta[index == 1]
    beta.logvar <- beta[index == 2]
    if (cate == 1){
      beta.gamma <- beta[index == 3]}
    if (cate == 2){
      beta.gamma <- rep(0,length(beta[index == 3]))}
    if (cate == 3){                       # change 2
      beta.mu <- rep(0,length(beta[index == 1]))
      beta.gamma <- beta[index == 3]}
    
    est <- array(NA, c(3,n,m))              
    #estimation of B %*% beta; NOT the true var and skewness
    est[1,,]  <- matrix(Basis.list[[1]] %*% beta.mu, n, m, byrow=F)
    est[2,,] <-  matrix(Basis.list[[2]] %*% beta.logvar, n, m, byrow=F)
    est[3,,] <-  matrix(c(Basis.list[[3]] %*% beta.gamma), n, m, byrow=T)               
    
    mu <- est[1,,]
    sigma <- exp(est[2,,]/2)
    x <- (dataone - mu)/sigma
    M <- 0.99527174643
    gamma1 <- (pnorm(est[3,,])*2 - 1)*M
    
    
    a2 <- D.lg(x,gamma1)                 #derivatives of log(g) wrt x and gamma1   
    g1 <- apply(c(-a2$D1/sigma) * Basis.list[[1]], 2, sum)         
    g2 <-  -1/2* apply(c(a2$D1*x +1)*Basis.list[[2]], 2, sum)   
    g3 <-   apply(a2$D2*2*M*dnorm(est[3,,]), 2, sum)             
    wtb.g3 <- apply((matrix(g3, nrow(Basis.list[[3]]), ncol(Basis.list[[3]])) * Basis.list[[3]]), 2, sum) 
    
    # calculate the negative gradient as gr
    if (cate ==1){  gr <- -c(g1,g2,wtb.g3)} 
    if (cate ==2){  gr <- -c(g1,g2,rep(0,length(beta[index == 3])))} 
    if (cate ==3){  gr <- -c(rep(0,length(beta[index == 1])),g2,wtb.g3)} 
    
    return(gr)}



### case2.b.initial: function to obtain the initial estimates of functional parameters (mean, var, shape and skewness) when case = 2.
# INPUT  
# @ y - observed data, in a form a matrix 
# @ timepoints, covariate for each subject (not in order); 
# @ nbasis.mean : number of bases when smoothing the mean;
# @ gam.method = smoothing method
# @ bin = the length of bin to estimate the variance 
# @ skew.method = estimation method for skewness, values = "mome" or "mle"
# @ cate = 1:all to be estimated; 2: shape fixed at 0 ; 3: mean to be fixed at 0.
# OUTPUT  
# a list of initial estimate of parmeters (length 4: mean, var, shape and skewness)          

case2.b.initial <- function(y, tp, cp, nbasis.mean = 10, gam.method = "REML", 
                            bin = 10, skew.method = "mle", 
                            cate = 1){   
  
  n.sub <- nrow(y)                            
  n.tim <- ncol(y)
  ######################################################## 
  ####initial estimates of mean via bivariate smooth #####
  ########################################################
  if (cate != 3){
    y.data <- as.vector(y)
    x1 <- rep(cp, n.tim)
    x2 <- rep(tp, each=n.sub)
    x.data <- data.frame(x1=x1,x2=x2,y.data=y.data)
    smooth.data <- gam(y.data~te(x1,x2, k=nbasis.mean), data=x.data, 
                       method = gam.method, 
                       omit.missing=TRUE )
    
    fitted <- predict(smooth.data,newdata=data.frame(x1=x1,x2=x2)) 
    
    #initial estimate for mean matrix by bivariate smoothing
    ini.mean <- matrix(fitted, ncol=n.tim)
  }  
  
  if (cate == 3){            
    ini.mean <- matrix(0, nrow=n.sub, ncol=n.tim)
  }
  
  ini.noise <- y - ini.mean   #initial recover of the random noise
  

  
  ########################################################## 
  ####initial estimates of variance via binning method #####      
  ##########################################################
  
  
      a <- ini.noise
      ini.var<- matrix(NA, nrow(a), ncol(a))
      for (i in 1:nrow(a)){
        north <- max(i - bin, 1)
        south <- min(i + bin, nrow(a))
        for (j in 1:ncol(a)){
          west <- max(j - bin,1)
          east <- min(j + bin, ncol(a) )           
          scope <- as.vector(a[north:south, west:east])
          ini.var[i,j] <- var(scope) }
      } 
                            
  
  ############################################################ 
  ####initial estimates of skewness via method of moment #####
  ############################################################   
       
  scale.data <- (y - ini.mean)/sqrt(ini.var) #scale data  
  if (cate != 2){
    if (skew.method == "mome"){
      # use sample skewness as initial estimates
      # require(moments)
      est.skew = as.vector(unlist(lapply(c(1:ncol(scale.data)), function(k) 
        skewness(scale.data[,k]))))   
      M <-  0.99527174643  
      # truncate the est of skewness by [-M,M]
      est.skew[abs(est.skew) >= M] <- 0.99 * sign(est.skew[abs(est.skew) >= M])  
    }
    
    if (skew.method == "mle"){
      est.skew = as.vector(unlist(lapply(1:n.tim, function(k) 
        coef(selm(scale.data[,k] ~ 1, family = "SN"))["gamma1"])))
    }
    
    est.shape <- sapply(est.skew, shape.dp)     
    ini.skew <- est.skew
    ini.shape <- est.shape
  }
  
  
  if (cate == 2){
    ini.skew <- rep(0,ncol(scale.data))
    ini.shape <-rep(0,ncol(scale.data)) 
  }       
  ini.list <- list(mean = ini.mean, var = ini.var, skew = ini.skew, shape = ini.shape)
  return(ini.list)}


                      




#######################################################################
######### transform cp parameters to beta's ############################
#######################################################################
## function cp2beta
## input:   cp.list - list of parameters with mean, var, skew, shape
##          Basis.list <- list of basis matrix with 3 components B.mean, B.logvar, B.gamma
## output:  a vector of beta's
cp2beta <- function(cp.list, Basis.list) {
  
  Basis.list$index <- NULL
  S.list <- lapply(Basis.list, function(B) ginv(t(B) %*% B) %*% t(B) ) #projection matrix 
  
  trans.cp.list <- as.list(c(1:3))  #transform cp.list to real line
  trans.cp.list[[1]] <- cp.list$mean
  trans.cp.list[[2]] <- log(cp.list$var)
  
  M <-  0.99527174643 
  trans.cp.list[[3]] <- qnorm((cp.list$skew/M + 1)/2)         #pnorm transformation
  
  
  beta.list <- lapply(c(1:3), function(k) S.list[[k]] %*% c(trans.cp.list[[k]]))
  #index <- lapply(c(1:3), function(k) rep(k,length(beta.list[[k]])))        
  return(unlist(beta.list))    
}


#######################################################################
######### beta2cp: transform beta's to cp parameters ##################
#######################################################################
## input:   vec.beta - vector of all the beta coefficients for mean, var, skew, shape
##          Basis.list <- list of basis matrix with 3 components B.mean, B.logvar, B.gamma
##                        will be modified to 3 components
## output:  a list of parameters (mean, var, skew, shape) 
beta2cp <- function(vec.beta, Basis.list) {
  index = unlist(lapply(c(1:3), function(k) rep(k, ncol(Basis.list[[k]]))))  
  beta.mu <- vec.beta[index == 1]
  beta.logvar <- vec.beta[index == 2]
  beta.gamma <- vec.beta[index == 3]
  
  gamma <- c(Basis.list[[3]] %*% beta.gamma)
  m <- length(gamma)
  
  est.mu <- matrix(Basis.list[[1]] %*% beta.mu, ncol=m, byrow=F)
  est.logvar <- matrix(Basis.list[[2]] %*% beta.logvar, ncol=m, byrow=F)
  est.gamma <- matrix(gamma, ncol=m, byrow=T)      
  est.var <- exp(est.logvar)   
  
  est.skew <- (pnorm(est.gamma)*2 - 1)*0.99527174643 
  est.shape <- sapply(est.skew, shape.dp)
  
  cp.list <- list(mean = est.mu, var=est.var, skew=est.skew, shape=est.shape)
  return(cp.list)}










# Function used for model estimation with given knots 

cSFM.est <- function(data, tp, cp, nknots.tp, nknots.cp, degree.poly = c(3, 3, 3), 
          method = c("cSFM", "cSFM0", "2cSFM"), bi.level = 2,  
          nbasis.mean = 10, gam.method = "REML"){
  
    n.sub <- nrow(data)
    n.tim <- ncol(data)
    my.data <- data        
    
    if (method == "cSFM") cate <- 1
    if (method == "cSFM0") cate <- 2
    if (method == "2cSFM") cate <- 3
    
    if (method == "2cSFM"){
      y.data <- as.vector(my.data)
    x1 <- rep(cp, n.tim)
    x2 <- rep(tp, each=n.sub)
    x.data <- data.frame(x1=x1,x2=x2,y.data=y.data)
    gam.mean <- gam(y.data~te(x1,x2, k=nbasis.mean), data=x.data, 
                       method = gam.method, 
                       omit.missing=TRUE )
    
    fitted <- predict(gam.mean,newdata=data.frame(x1=x1,x2=x2)) 
    
     # estimates for mean matrix by bivariate smoothing
    est.mean <- matrix(fitted, ncol=n.tim)
      # use residuals as the new input
    my.data <- my.data - est.mean
    }
    
    # obtain initial estimates for the three parameters. 
    if (bi.level == 2){
   ini.orig2 <- case2.b.initial(my.data, tp, cp, cate=cate)
   }

#     if (bi.level == 0){
#    ini.orig2 <- case1.initial(my.data)}
    
    # generate the basis system for all three paraemter functions
    if (bi.level == 2){
      sub.case = c(2,2,1)
    }
    B.mu <- kpbb(tp, cp, nknots.tp = nknots.tp[1], 
                       nknots.cp = nknots.cp[1], 
                       sub.case = sub.case[1], degree.poly = degree.poly[1])
    B.logvar <- kpbb(tp, cp, nknots.tp = nknots.tp[2], 
                           nknots.cp = nknots.cp[2],
                           sub.case = sub.case[2], degree.poly = degree.poly[2])
    B.gamma <- kpbb(tp, cp, nknots.tp = nknots.tp[3], 
                          nknots.cp = nknots.cp[3], 
                          sub.case = sub.case[3], degree.poly = degree.poly[3])    
    sm.basis = list(sm.mu = B.mu, sm.logvar=B.logvar, sm.gamma = B.gamma)
     
    # initial coefficients
    ini.beta2 <- cp2beta(ini.orig2, sm.basis)  
    
    # optimal coefficients
    beta.hat <- optim(par=ini.beta2, fn = case2.unmll.optim, gr=case2.gr, 
                        dataone=my.data, Basis.list=sm.basis, cate=cate, 
                        method="BFGS",control=list(trace=0,  maxit=5000)
                        )
    
    if (beta.hat$convergence != 0 ) print(beta.hat$message)
    cp.est <- beta2cp(beta.hat$par, sm.basis)     
    
    if (method == "2cSFM"){
      cp.est$mean <- est.mean
    }
    
    # obtain copula and AIC
    copula <- matrix(NA, nrow(data), ncol(data))
    for (l1 in 1:nrow(data)){
      for (l2 in 1:ncol(data)){
        temp.dp <- cp2dp(c(cp.est$mean[l1,l2], 
                              sqrt(cp.est$var[l1,l2]), cp.est$skew[l2]), family = "SN")
        copula[l1,l2] <- psn(data[l1,l2], dp = temp.dp)
      }
    }
    back.z <- qnorm(copula)   # latent Gaussian process
    x <- (data - cp.est$mean)/sqrt(cp.est$var)    
    tau <- sin(cor(x,method="kendall")*pi/2)    
    est.corr <- tau    
    lstar <- rep(NA, nrow(data))    
    for (l3 in 1:nrow(data)){
      data.z <- matrix((back.z[l3,] ), ncol(data),1)
      lstar[l3] <- -0.5*ncol(data)*log(2*pi) - 0.5*det(est.corr) - as.vector(0.5 * t(data.z) %*% solve(est.corr) %*% (data.z)) }        
    
    parlength <-  length(beta.hat$par)    
    AIC = 2*beta.hat$value + 2*parlength -2*sum(lstar) #AIC with dependence     
    
    # save estimation objects and return
    cp.hat <- cp.est
    corr.est <- tau
    num.pars <- parlength
    logL1 <- -beta.hat$value
    logL2 <- sum(lstar)
    ret.objects = c("beta.hat", "cp.hat", "copula", "corr.est", "AIC", 
                    "num.pars", "logL1", "logL2", "sm.basis")
    
    if (method == "2cSFM"){
      ret.objects <- c(ret.objects, "gam.mean")
    }
    
   ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
   names(ret) = ret.objects
    
   ret$call <- match.call()
    
   class(ret) = "cSFM"
   attributes(ret)$method <- method
   return(ret)    
}

# functions for print
print.cSFM <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients (mean):\n")
  print(x$beta.hat$par[x$sm.basis$index == 1])
  cat("\nCoefficients (logvar):\n")
  print(x$beta.hat$par[x$sm.basis$index == 2])
  cat("\nCoefficients (gamma):\n")
  cat("\nAIC = ", x$AIC)
}

# functions to extract fitted values
fitted.cSFM <- function(object, quantile = TRUE, quantile.level = c(.50, .80, .90, .95, .99), ... )
{ 
  x <- object
  est.col <- with(x, 
                  list(mean = cp.hat$mean, logvar = log(cp.hat$var), skew = cp.hat$skew)
                  )
  
  # extract quantile estimates
  
  
  if (quantile == TRUE){
    est.quantile = array(0, dim = c(dim(est.col$mean), length(quantile.level)))  
    for(p in 1:length(quantile.level)){
      
      for (j in 1:ncol(est.col$mean)){
        temp2 <- qsn(quantile.level[p], dp = cp2dp(c(0,1,est.col$skew[j]), family = "SN"))
        est.quantile[,j,p] <- est.col$mean[,j] + exp(est.col$logvar[,j]/2) * temp2      
      }
    }
  }
  est.col$quantile <- est.quantile
  est.col
}


# functions for Predication 


predict.kpbb <- function(object, tp.valid, cp.valid = NULL, ...)
{
  sub.case = ("tp" %in% names(attributes(object))) + 1
  
  if (sub.case == 1)
  {
    B.tim <- predict(object, tp.valid) 
    return(B.tim)
  }
  
  if (sub.case == 2)
  {
    a <- c(list(x = tp.valid), attributes(object)$tp[c("degree", "knots",
                                                       "Boundary.knots", "intercept")])
    B.tim <- do.call("bs", a)
    
    b <- c(list(x = cp.valid), attributes(object)$cp[c("degree", "knots",
                                                       "Boundary.knots", "intercept")])
    B.cov <- do.call("bs", b)
    
    B.tensor <- kronecker(B.tim, B.cov)
    return(B.tensor)    
    
  }
  
}




## regular FPCA for univariate with given covariance matrix or sample covariance matirix; the curves could be centered or assumbed to be centerred
## Input: 
# @ Y: fully observed curves Y
# @ Kendall: given covariance matrix for Y (raw matrix up to smoothing)
# @ Y.pred: partically observed curves
# @ center: if true, then we center the curves; otherwise, assume the curves to be centered already

uni.fpca <- function (Y, Kendall=NULL, Y.pred= NULL, nbasis.mean = 10, nbasis = 10,
                      pve = 0.99, npc = NULL, center = TRUE, gam.method = "REML"){
  
  # require(mgcv)
  # require(MASS)
  if (is.null(Y.pred)) 
    Y.pred = Y
  D = NCOL(Y)
  I = NROW(Y)
  I.pred = NROW(Y.pred)
  d.vec = rep(1:D, each = I)
  
  if (center){
    gam0 = gam(as.vector(Y) ~ s(d.vec, k = nbasis.mean), method=gam.method)
    mu = predict(gam0, newdata = data.frame(d.vec = 1:D))
    Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)}
  
  if (!center){
    mu = rep(0, D)
    Y.tilde = Y
  }
  
  if (is.null(Kendall)){    
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    for (i in 1:I) {
      obs.points = which(!is.na(Y[i, ]))
      cov.count[obs.points, obs.points] = cov.count[obs.points, 
                                                    obs.points] + 1
      cov.sum[obs.points, obs.points] = cov.sum[obs.points, 
                                                obs.points] + 
                                                  tcrossprod(Y.tilde[i, obs.points])
    }
    G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)}
  
  if (!is.null(Kendall)) G.0 = Kendall
  
  diag.G0 = diag(G.0)
  diag(G.0) = NA
  
  row.vec = rep(1:D, each = D)
  col.vec = rep(1:D, D)
  npc.0 = matrix(predict(gam(as.vector(G.0) ~ te(row.vec, 
                                                 col.vec, k = nbasis), method=gam.method), 
                         newdata = data.frame(row.vec = row.vec, 
                                              col.vec = col.vec)), D, D)
  npc.0 = (npc.0 + t(npc.0))/2
  
  
  
  evalues = eigen(npc.0, symmetric = TRUE, only.values = TRUE)$values
  evalues = replace(evalues, which(evalues <= 0), 0)
  npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > 
    pve)), npc)
  efunctions = matrix(eigen(npc.0, symmetric = TRUE)$vectors[, 
                                                             seq(len = npc)], nrow = D, ncol = npc)
  evalues = eigen(npc.0, symmetric = TRUE, only.values = TRUE)$values[1:npc]
  
  cov.hat = efunctions %*% tcrossprod(diag(evalues, nrow = npc, 
                                           ncol = npc), efunctions)
  K.tilde = cov.hat
  DIAG = (diag.G0 - diag(cov.hat))[floor(D * 0.2):ceiling(D * 
    0.8)]
  sigma2 = max(mean(DIAG, na.rm = TRUE), 0)
  K.hat = K.tilde
  diag(K.hat) = diag(K.hat) + sigma2
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)
  Z = efunctions
  Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
  Yhat = matrix(0, nrow = I.pred, ncol = D)
  scores = matrix(NA, nrow = I.pred, ncol = npc)    
  for (i.subj in 1:I.pred) {
    obs.points = which(!is.na(Y.pred[i.subj, ]))
    if (sigma2 == 0 & length(obs.points) < npc) {
      stop("Measurement error estimated to be zero 
                 and there are fewer observed points than PCs; 
           scores cannot be estimated.")
        }
    Zcur = matrix(Z[obs.points, ], nrow = length(obs.points), 
                  ncol = dim(Z)[2])
    ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)
    scores[i.subj, ] = ZtZ_sD.inv %*% t(Zcur) %*% (Y.tilde[i.subj, 
                                                           obs.points])
    Yhat[i.subj, ] = t(as.matrix(mu)) + scores[i.subj, ] %*% 
      t(efunctions)
    
  }
  ret.objects = c("Yhat", "scores", "mu", "efunctions", "evalues", 
                  "npc", "K.tilde", "K.hat", "sigma2")
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  return(ret)
}



## predication 
## create a function with input:
# object: an object of estimates returned by the function est.cSFM with class "HAC"

# @ newdata: patically observed curves matrix to be predicted
# @ cp.valid : covariate for Y.p
# @ tp.valid : timepoints for all curves; shared by Y.f and Y.p

predict.cSFM <- function (object, newdata, cp.valid, tp.valid = NULL, ...) {
  
  
  y.p <- newdata
  if (is.null(tp.valid))
  {
    n.time = ncol(object$copula)
    tp.valid = seq(from = 0, to = 1, length = n.time)
  }
  
  
  # generate the basis for predication
  B.mu <- predict(object$sm.basis$sm.mu, tp.valid, cp.valid)
  B.logvar <- predict(object$sm.basis$sm.logvar, tp.valid, cp.valid)
  B.gamma <- predict(object$sm.basis$sm.gamma, tp.valid, cp.valid)
  index <- c(rep(1, ncol(B.mu)), rep(2, ncol(B.logvar)), rep(3, ncol(B.gamma)))
  bi.level <- object$sm.basis$bi.level
  sm.basis <- list(sm.mu = B.mu, sm.logvar=B.logvar, sm.gamma = B.gamma, 
                   index = index, bi.level = bi.level)
  
  cp.p <- beta2cp(object$beta.hat$par, sm.basis)
  
  if (attr(object, "method") == "2cSFM"){
    # add back the bivariate mean for 2cSFM method 
    #predicated mean matrix by bivariate smoothing
    mu.p <- predict(object$gam.mean,
                    newdata=data.frame(x1=rep(cp.valid, length(tp.valid)),
                                       x2=rep(tp.valid, each=length(cp.valid))))  
    mu.p <- matrix(mu.p, ncol = length(tp.valid))  
    cp.p$mean <- mu.p # add back the mean estimates
  }    
  
  ## recover copula   
  scale.y.p = (y.p - cp.p$mean)/sqrt(cp.p$var)
  copula.p = matrix(NA, nrow(scale.y.p), ncol(scale.y.p))
    for (ith in 1:nrow(scale.y.p)){
      copula.p[ith, ] <- sapply(c(1:ncol(scale.y.p)), 
                                function(k) ifelse(is.na(scale.y.p[ith, k]), NA, 
                                                   psn(scale.y.p[ith, k], 
                                                       dp = cp2dp(c(0, 1, cp.p$skew[k]), family = "SN"))))
    }
  
  GP.p = matrix(sapply(copula.p, function(k) ifelse(is.na(k), NA, qnorm(k))), nrow(copula.p), ncol(copula.p))
  
  # robust mean and sd for GP.p
  rob.mu.sigma = apply(GP.p, 1, function(k) hubers(k[is.finite(k)]))
  rob.mu = sapply(rob.mu.sigma, function(k) k$mu)
  rob.sigma = sapply(rob.mu.sigma, function(k) k$s)
  
  bound = 5 # may try 4 and 6 # report how many are extreme 
  up.bound = unlist(lapply(rob.mu.sigma, function(k) k$mu + bound*k$s))
  low.bound = unlist(lapply(rob.mu.sigma, function(k) k$mu - bound*k$s))
  
  truc.GP <- function(GP.vector, lower, upper){
    GP.vector[GP.vector< lower | GP.vector == -Inf] = lower
    GP.vector[GP.vector > upper | GP.vector == Inf] = upper
    GP.vector
  }
  
  t.GP.p = sapply(c(1:ncol(GP.p)), function(k) truc.GP(GP.p[,k], low.bound[k], up.bound[k]))
  
  
  # use kendall's tau to smooth t.GP.p
  Y = qnorm(object$copula)  
  Kendall =  object$corr.est 
  Y.pred = t.GP.p # use truncated GP.p 
  
  # Functional PCA given the covariance matrix as Kendall's Tau  
  k = uni.fpca(Y, Kendall, Y.pred)  
  # adjust the predicated GP  when recoving copula  
  std.z = (k$Yhat - matrix(k$mu, nrow=nrow(k$Yhat), ncol=ncol(k$Yhat), byrow=T))/ matrix(sqrt(diag(k$K.hat)), nrow=nrow(k$Yhat),ncol=ncol(k$Yhat), byrow=T)
  copula.pred = pnorm(std.z)
  
  # recover the latent process from corpula
  sn.process <- matrix(NA, nrow(copula.pred), ncol(copula.pred))
  for (ith in 1:nrow(copula.pred)){
    sn.process[ith,] <- sapply(c(1:ncol(copula.pred)), 
                               function(k) qsn(copula.pred[ith, k],
                                               dp=cp2dp(c(0,1,cp.p$skew[k]), family = "SN"))) 
  }
  
  # prediction matrix
  yhat.pred = cp.p$mean + sqrt(cp.p$var) * sn.process
  
  return(yhat.pred)
}


# knots selection 

cSFM.est.parallel <- function(data, tp, cp, nknots.tp = NULL, nknots.cp = NULL, max.knots = NULL, 
                              parallel = TRUE, num.core = NULL, method = c("cSFM", "cSFM0", "2cSFM"), bi.level = 2,
                              degree.poly = c(3, 3, 3), nbasis.mean = 10, gam.method = "REML"){
  
  if (!is.null(max.knots)){
    a <- expand.grid (nk1 = 1:max.knots, nk2 = 1:max.knots, nk3 = 1:max.knots)
    index1 <- with(a, (nk1 > nk2) & (nk2 > nk3))
    all.knots <- a[index1,]
    knots.mat <- as.matrix(all.knots)
    knots.mat <- cbind(knots.mat, knots.mat)   
  } else {
    knots.mat = cbind(nknots.tp, nknots.cp)
  }
  
  colnames(knots.mat) <- c("nk1.tp", "nk2.tp", "nk3.tp", "nk1.cp", "nk2.cp", "nk3.cp")
  rownames(knots.mat) <- c(1:nrow(knots.mat))
  
  cat(nrow(knots.mat), "knots settings are considerred. \n")
  temp.fun <- function(i) {
    cSFM.est(data, tp, cp, knots.mat[i, 1:3], knots.mat[i, 4:ncol(knots.mat)], degree.poly = degree.poly, 
             method = method, bi.level = bi.level, 
             nbasis.mean = nbasis.mean, gam.method = gam.method)
  }
  
  start <- proc.time()
  cat("Knots selection begins ... ")
  if (!parallel){
    cat("parallel computing not used ...")
    temp <- lapply(as.list(1:nrow(knots.mat)), temp.fun)
  }
  
  if (parallel){
    # require(parallel)
    cat("parallel computing applied ...")
    if (is.null(num.core)){
      num.core <- detectCores()
    }
    
    cat(num.core, "Cores are used ... ")
    cl <- makeCluster (getOption ("cl.cores", num.core))
    
    clusterEvalQ(cl, require(cSFM)) 
    clusterExport (cl, c("data", "tp", "cp", "knots.mat", "degree.poly", "method", "bi.level", "nbasis.mean", "gam.method"), 
                   envir = environment())
    temp <- clusterMap (cl, temp.fun, 1:nrow(knots.mat), .scheduling = "dynamic")
    stopCluster(cl)
  }
  duration = proc.time() - start
  cat("Done! Took", duration[3], "s \n")
  
  
  AIC.vector <- unlist(lapply(temp, function(k) k$AIC))
  
  best.index = which.min(AIC.vector)
  ret = list(best.cSFM = temp[[best.index]], knots.mat = knots.mat, AIC = AIC.vector)
  return(ret)
  
}
