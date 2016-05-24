BTLLasso.ctrl <- 
function (adaptive = TRUE, norm = c("L1","L2"), epsilon = 1e-4, lambda2 = 1e-4,
          c = 1e-9, penal.diffs = TRUE, return.design = TRUE) 
{
  norm <- match.arg(norm)
  RET <- list(adaptive = adaptive, norm = norm, epsilon = epsilon, lambda2 = lambda2,
             c = c, penal.diffs = penal.diffs, return.design = return.design)
  RET
}


BTLLasso <- function(Y, X, lambda, control = BTLLasso.ctrl(), trace = TRUE){
  
  vardiffs <- abs(apply(X,2,var)-1)
  
  
  if(!is.matrix(X))
    stop("X has to be a matrix")
  
  if(!is.matrix(Y))
    stop("Y has to be a matrix")

  n <- nrow(Y)
  I <- ncol(Y)
  m <- (1 + sqrt(1+8*I))/2
  p <- ncol(X)
  
  
  ### extract all control arguments
  adaptive <- control$adaptive
  norm <- control$norm
  epsilon <- control$epsilon
  lambda2 <- control$lambda2
  c <- control$c
  penal.diffs <- control$penal.diffs
  return.design <- control$return.design
  ###
  
  labels.raw <- colnames(Y)[1:(m-1)]
  labels <- c(word(labels.raw)[1],word(labels.raw,3))
  labels <- labels[c(length(labels),1:(length(labels)-1))]
  
  Y <-as.ordered(c(Y))
  
  q <- length(levels(Y))-1
  k <- q+1
  get.resp <- function(x){
    as.numeric(as.numeric(x) <= 1:q)
  }
  

  
  index <- rep(1:p,each=m-1)

  
  acoefs <- diag(p*(m-1))
  
  if(penal.diffs){
    help.pen <- matrix(0,ncol=choose(m-1,2),nrow=m-1)
    combis <- combn(m-1,2)
    for(ff in 1:ncol(combis)){
      help.pen[combis[1,ff],ff] <- 1
      help.pen[combis[2,ff],ff] <- -1
    }
    for(pp in 1:p){
      m.above <- matrix(rep(matrix(0,ncol=choose(m-1,2),nrow=m-1),pp-1),ncol=choose(m-1,2))
      m.below <- matrix(rep(matrix(0,ncol=choose(m-1,2),nrow=m-1),p-pp),ncol=choose(m-1,2))
      acoefs <- cbind(acoefs,rbind(m.above,help.pen,m.below))
    } 
  }
  
  
  acoefs <- rbind(matrix(0,nrow=floor(q/2)+m-1,ncol=ncol(acoefs)),acoefs)
  

  
  response <- c()
  for(i in 1:length(Y)){
    response <- c(response,get.resp(Y[i]))
  }

  
  design2 <- create.design(X,m)
  design <- t(matrix(rep(c(design2),each=q),nrow=ncol(design2),byrow=TRUE))
  
  n.theta <- floor(q/2)
  
  if(k>2){
  theta.design <- matrix(0,ncol=n.theta,nrow=nrow(design))
  
  for(i in 1:n.theta){
    vec1 <- rep(0,q)
    vec1[c(i,q-i+1)] <- c(1,-1)
    theta.design[,i] <- rep(vec1,nrow(design2))
  }
  
  design <- cbind(theta.design, design)
  }


  na.response <- which(is.na(response)) 
  na.design <- which(apply(design,1,function(x){any(is.na(x))}))
  na.index <- unique(c(na.response, na.design))
  na.Y <- (na.index/q)[na.index%%q == 0]

  if(length(na.index)>0){
  response <- response[-na.index]
  design <- design[-na.index,]
  Y <- Y[-na.Y]
  }
  

  coefs <- matrix(0,nrow=length(lambda),ncol=ncol(design))
  df <- c()
  start <- NULL
  
  if(adaptive){
    if(k>2){
    m0 <- cum.fit.Cpp(response, design, kat=k, epsilon = epsilon, start=start, 
                  acoefs=acoefs, lambda=lambda2, max.iter=100, norm = norm,
                  adaptive = NULL, control = list(c = c, gama = 20, index = index), 
                  m = m, hat.matrix = FALSE, lambda2 = lambda2)
    }else{
      m0 <- bin.fit.Cpp(response, design, kat=k, epsilon = epsilon, start=start, 
                        acoefs=acoefs, lambda=lambda2, max.iter=100, norm = norm,
                        adaptive = NULL, control = list(c = c, gama= 20, index = index), 
                        m = m, hat.matrix = FALSE, lambda2 = lambda2)
    }
    adaptive <- m0$coef
  }else{
    adaptive <- NULL
  }
  
  for(i in seq_along(lambda)){
    if(trace){cat("lambda =",lambda[i],"\n")}
    if(k>2){
  m1 <- cum.fit.Cpp(response, design, kat=k, epsilon = epsilon, start=start, 
                  acoefs=acoefs, lambda=lambda[i], max.iter=100, norm = norm,
                  adaptive = adaptive, control = list(c = c, gama = 20, index = index), 
                  m = m, hat.matrix = FALSE, lambda2 = lambda2)
    }else{
      m1 <- bin.fit.Cpp(response, design, kat=k, epsilon = epsilon, start=start, 
                        acoefs=acoefs, lambda=lambda[i], max.iter=100, norm = norm,
                        adaptive = adaptive, control = list(c = c, gama = 20, index = index), 
                        m = m, hat.matrix = FALSE, lambda2 = lambda2)
    }
  coefs[i,] <- m1$coef
  start <- m1$coef
    df[i] <- m1$df
 
  }
  
  
  logLik <- c()
  for(j in 1:nrow(coefs)){
    logLik[j] <- loglik(coefs[j,],Y,design,k)
  }

  
  gamma <- coefs[,-(1:(n.theta+m-1))]
  
  if(norm=="grouped"){
  norm.gamma <- matrix(0, ncol = nrow(gamma), nrow = p)
  for (j in 1:nrow(gamma)) {
    gamma.lambda <- matrix(gamma[j,],ncol=m-1,byrow=TRUE)
    norm.gamma[,j] <- sqrt(rowSums(gamma.lambda^2))
  }

  norm.gamma <- norm.gamma/norm.gamma[,ncol(norm.gamma)]
  
  df2 <- n.theta + m-1 + colSums(norm.gamma!=0) + colSums(norm.gamma)*(m-2)
  }else{
    df2 <- n.theta + m-1 + rowSums(gamma!=0)
  }
  
  
  aic <- -2*logLik + 2 *df
  bic <- -2*logLik + log(nrow(design)) *df
  
  aic2 <- -2*logLik + 2 *df2
  bic2 <- -2*logLik + log(nrow(design)) *df2
  
  if(!return.design){design <- NULL}
  
  coefs.repar <- round(expand.coefs(coefs,m=m,n.theta=n.theta),3)
  
  ret.list <- list(coefs = coefs, coefs.repar = coefs.repar, logLik = logLik, 
                   design = design, Y = Y, q = q, acoefs = acoefs, response = response, 
                   n = n, I = I, m = m, p = p, X = X, n.theta = n.theta, lambda = lambda, 
                   labels = labels, epsilon = epsilon)

  class(ret.list) <- "BTLLasso"
  
  return(ret.list)
  
}