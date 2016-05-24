expit <- function(x){ return( exp(x)/(exp(x)+1) ) }

#' @name datagenerator
#' @title simulate data
#' @description Simulate data for demonstration of flarcc.
#' @param n sample size for each study, a vector of length \code{K}, the number of studies; can also be an scalar, to specify equal sample size
#' @param beta0 coefficient matrix, with dimension \code{K * p}, where K is the number of studies and p is the number of covariates 
#' @param family "gaussian" for continuous response, "binomial" for binary response, "poisson" for count response
#' @param seed set random seed for data generation
#' @keywords data generator
#' @details These data sets are artifical, and used to test out some features of flarcc.
#' @return a simulated data frame will be returned
#' @export
# data generator, with covariates from N(0,1)
datagenerator <- function(n=n, beta0=beta0, family="gaussian", seed=seed){
  set.seed(seed+1111)
  p <- ncol(beta0); K <- nrow(beta0); 
  if(length(n)!=1 & length(n)!=K){ stop("n has to be scalar or a vector of length K = nrow(beta0).") }
  groupvec <- rep(c(1:K), times=n); N <- length(groupvec);
  betamat <- matrix(NA, N, p); for(i in 1:p){ betamat[,i] <- rep(beta0[,i], times=n) }
  if(p==2){ X <- matrix(c(rep(1,N), rnorm((p-1)*N)), N, p) }
  if(p>=3){ 
    sigma <- matrix(0.3, p-1, p-1); diag(sigma) <- 1
    X <- cbind(rep(1,N), mvrnorm(N, mu=rep(0, (p-1)), Sigma=sigma))
  }
  if(family=="gaussian"){ y <- rowSums(X*betamat) + rnorm(N) }
  if(family=="binomial"){ y <- rbinom(N, 1, expit(rowSums(X*betamat))) }
  if(family=="poisson"){  y <- rpois(N, exp(rowSums(X*betamat))) }
  colnames(X) <- paste("x", c(0:(p-1)), sep="")
  data <- data.frame(y=y,X,group=groupvec); data <- data[,-2] 
  data <- data[sample(1:N, N, replace=F),]; row.names(data) <- NULL
  rm(betamat, i, n, p, K, N, groupvec, seed, y, family, beta0)
  return(data)
}

#' @name metafuse
#' @title fit a GLM with fusion penalty for data integraion
#' @description Fit a GLM with fusion penalty on coefficients within each covariate, generate solution path for model selection.
#' @param X a matrix (or vector) of predictor(s), with dimensions of N*p, where N is the total sample size of all studies
#' @param y a vector of response, with length N, the total sample size of all studies
#' @param sid study id, numbered from 1 to K
#' @param fuse.which a vector of a subset of integers from 0 to p, indicating which covariates to be considered for fusion; 0 corresponds to intercept
#' @param family "gaussian" for continuous response, "binomial" for binary response, "poisson" for count response
#' @param intercept if TRUE, intercept will be included in the model
#' @param alpha the ratio of sparsity penalty to fusion penalty, default is 0 (no penalty on sparsity)
#' @param criterion "AIC" for AIC, "BIC" for BIC, "EBIC" for extended BIC
#' @param verbose if TRUE, output fusion events and tuning parameter lambda 
#' @param plots if TRUE, create plots of solution paths and clustering trees
#' @param loglambda if TRUE, lambda will be plot in log-10 transformed scale
#' @details More details to include.
#' @return a data matrix containing coefficient estimates, cluster size, optimal lambda by BIC/EBIC, and friction of fusion for each covariate 
#' @export
metafuse <- function(X=X, y=y, sid=sid, fuse.which=c(0:ncol(X)), family="gaussian", intercept=TRUE, alpha=0,
                     criterion="EBIC", verbose=TRUE, plots=TRUE, loglambda=TRUE){
  # argument check
  if(missing(X)){ stop("Need to specify X.") }
  if(missing(y)){ stop("Need to specify y.") }
  if(missing(sid)){ stop("Need to specify sid, the study ID for each subject") }
  if(min(fuse.which)<0 | max(fuse.which)>ncol(X)){ stop("Argument fuse.which must be an ordered ingeter vector from 0 to ncol(X).") }
  if(intercept==F & any(fuse.which==0)){ stop("No intercept to fuse, argument fuse.which cannot contain 0.") }
  
  # only complete cases
  complete <- complete.cases(cbind(X,y,sid))
  X        <- subset(X,   subset=complete)
  y        <- subset(y,   subset=complete)
  sid      <- subset(sid, subset=complete)
  
  studies  <- sort(unique(sid))
  K        <- length(studies)
  n.k      <- sapply(studies, function(x){sum(sid==x)})
  N        <- length(y)
  if(intercept==T){
    if(is.null(X)){X <- data.frame(Intercept=rep(1,N))
    } else        {X <- data.frame(Intercept=rep(1,N), X)}
  }
  X.names    <- names(X)
  p          <- ncol(X)
  fuse.which <- fuse.which+1*intercept
  weights    <- sapply(sid, function(x){((N/K) / n.k)[x]})
  
  X.b      <- X_trans(X=X, fuse.which=fuse.which, K=K, X.names=X.names, sid=sid, studies=studies)
  
  beta.glm <- glm_wrap(X=X.b, y=y, family=family, weights=weights)
  
  parg     <- pargroup(p=p, fuse.which=fuse.which, K=K)
  S        <- Sgen(p=p, parg=parg, fuse.which=fuse.which, K=K, beta=beta.glm)
  B        <- Bgen(p=p, parg=parg, fuse.which=fuse.which, K=K, beta=S%*%beta.glm)
  
  theta.glm<- drop(B %*% S %*% beta.glm)
  X.t      <- X.b %*% solve(B%*%S)
  
  center   <- centerlist(p=p, parg=parg)
  al.r     <- 1
  al.weights <- 1/(abs(theta.glm)^al.r)
  al.weights[center] <- alpha * al.weights[center]
  
  currentDF <- betaDF(beta=beta.glm, parg=parg, p=p, X.names=X.names)
  if(verbose==T){ 
    cat("************************** verbose (start) **************************\n")
    cat(paste("Lambda = ", round(0.00,4), "\n", sep="")); print(currentDF) 
  }
  groupPath <- getGrouping(beta.glm, parg=parg, p=p, fuse.which=fuse.which)
  lambdaCP  <- c(0)
  
  # list of summaries
  solPath <- c(); lamlist <- c(); dflist <- c()
  
  l <- 0; ll <- -4; dd <- 0.05
  while(!all(currentDF == 1)){
    
    if(all(X.t[,1]==1)){
      glmnet.ada.fit <- glmnet(X.t[,-1], y, family=family, lambda=l, weights=weights, intercept=T, standardize=F, penalty.factor=al.weights[-1])
      thetahat       <- drop(coef(glmnet.ada.fit))
    } else{
      glmnet.ada.fit <- glmnet(X.t, y, family=family, lambda=l, weights=weights, intercept=F, standardize=F, penalty.factor=al.weights)
      thetahat       <- drop(coef(glmnet.ada.fit)[-1])
    }
    betahat   <- drop(solve(B%*%S)%*%thetahat)
    
    tempDF <- betaDF(beta=betahat, parg=parg, p=p, X.names=X.names)
    if(any(currentDF != tempDF)) {
      currentDF <- tempDF
      if(verbose==T){ cat(paste("Lambda = ", round(l,4), "\n", sep="")); print(currentDF) }
      groupPath <- cbind(groupPath, getGrouping(betahat, parg=parg, p=p, fuse.which=fuse.which))
      lambdaCP  <- c(lambdaCP, l)
    }
    solPath   <- cbind(solPath, as.vector(betahat))
    lamlist   <- c(lamlist,l)
    
    l        <- 10^(ll)
    ll       <- ll+dd
    dflist   <- cbind(dflist, tempDF)
  }
  if(verbose==T){ cat("************************** verbose (end) **************************\n") }
  
  critPath          <- apply(solPath, 2, critfunc, family=family, criterion=criterion, N=N,X.b=X.b,y=y,K=K,p=p,sid=sid,parg=parg,fuse.which=fuse.which)
  colnames(solPath) <- round(lamlist,5); rownames(solPath) <- names(beta.glm)
  colnames(dflist)  <- round(lamlist,5)
  names(critPath)   <- round(lamlist,5)
  
  lambda_opt    <- lamlist[which.min(critPath)]
  lambda_opt_t  <- log(lambda_opt+1, base=10)
  
  beta_opt      <- solPath[,which.min(critPath)]
  betaopt_out   <- betahatReformat(beta=beta_opt, parg=parg, p=p, K=K)
  
  beta_opt_DF   <- dflist[,which.min(critPath)]
  
  lambda_fuse   <- lamlist[apply(dflist, 1, function(x){which(x==1)[1]})]
  lambda_fuse_t <- log(lambda_fuse+1, base=10)
  
  friction      <- sapply(lambda_fuse, function(x){ 1 - min(x, lambda_opt)/x })
  friction_t    <- sapply(lambda_fuse_t, function(x){ 1 - min(x, lambda_opt_t)/x })
  
  if.fuse       <- rep(0,p); if.fuse[fuse.which] <- 1; names(if.fuse) <- X.names
  
  out.beta      <- betaopt_out; colnames(out.beta) <- X.names
  out.info.o    <- rbind(DF=beta_opt_DF, lambda_opt, lambda_fuse, friction); colnames(out.info.o) <- X.names
  out.info.t    <- rbind(DF=beta_opt_DF, lambda_opt=lambda_opt_t, lambda_fuse=lambda_fuse_t, friction=friction_t); colnames(out.info.t) <- X.names
  out.info      <- list(original_scale=out.info.o, log10_scale=out.info.t)
  
  out           <- list(family=family, criterion=criterion, alpha=alpha, if.fuse=if.fuse, betahat=out.beta, betainfo=out.info)
  
  # generate all the plots
  if(plots==T){
    plotSolPath(solPath=solPath, lamlist=lamlist, parg=parg, p=p, lambda_opt=lambda_opt, loglambda=loglambda, X.names=X.names, fuse.which=fuse.which)
    readkey()
    plotCrit(critPath=critPath, lamlist=lamlist, lambda_opt=lambda_opt, loglambda=loglambda, criterion=criterion)
    for(plot.which in fuse.which){
      readkey()
      fusiogram(groupPath=groupPath[parg==plot.which,], lambdaCP=lambdaCP, beta=beta_opt, K=K, parg=parg,
                plot.which=plot.which, lambda_opt=min(lambda_opt,lambda_fuse[plot.which]),
                loglambda=loglambda, X.names=X.names)
    }
  }
  
  return(out)
}


#' @name metafuse.l
#' @title fit a GLM with fusion penalty for data integraion
#' @description Fit a GLM with fusion penalty on coefficients within each covariate at given lambda.
#' @param X a matrix (or vector) of predictor(s), with dimensions of N*p, where N is the total sample size of all studies
#' @param y a vector of response, with length N, the total sample size of all studies
#' @param sid study id, numbered from 1 to K
#' @param fuse.which a vector of a subset of integers from 0 to p, indicating which covariates to be considered for fusion; 0 corresponds to intercept
#' @param family "gaussian" for continuous response, "binomial" for binary response, "poisson" for count response
#' @param intercept if TRUE, intercept will be included in the model
#' @param alpha the ratio of sparsity penalty to fusion penalty, default is 0 (no penalty on sparsity)
#' @param lambda tuning parameter for fusion penalty
#' @details More details to include.
#' @return a data matrix containing coefficient estimates, cluster size, optimal lambda by BIC/EBIC, and friction of fusion for each covariate 
#' @export
metafuse.l <- function(X=X, y=y, sid=sid, fuse.which=c(0:ncol(X)), family="gaussian", intercept=TRUE, 
                       alpha=0, lambda=lambda){
  # argument check
  if(missing(X)){ stop("Need to specify X.") }
  if(missing(y)){ stop("Need to specify y.") }
  if(missing(sid)){ stop("Need to specify sid.") }
  if(min(fuse.which)<0 | max(fuse.which)>ncol(X)){ stop("Argument fuse.which must be an ordered ingeter vector from 0 to ncol(X).") }
  if(intercept==F & any(fuse.which==0)){ stop("No intercept to fuse, argument fuse.which cannot contain 0.") }
  
  # only complete cases
  complete <- complete.cases(cbind(X,y,sid))
  X        <- subset(X,   subset=complete)
  y        <- subset(y,   subset=complete)
  sid      <- subset(sid, subset=complete)
  
  studies  <- sort(unique(sid))
  K        <- length(studies)
  n.k      <- sapply(studies, function(x){sum(sid==x)})
  N        <- length(y)
  if(intercept==T){
    if(is.null(X)){X <- data.frame(Intercept=rep(1,N))
    } else        {X <- data.frame(Intercept=rep(1,N), X)}
  }
  X.names    <- names(X)
  p          <- ncol(X)
  fuse.which <- fuse.which+1*intercept
  weights    <- sapply(sid, function(x){((N/K) / n.k)[x]})
  
  X.b      <- X_trans(X=X, fuse.which=fuse.which, K=K, X.names=X.names, sid=sid, studies=studies)
  
  beta.glm <- glm_wrap(X=X.b, y=y, family=family, weights=weights)
  
  parg     <- pargroup(p=p, fuse.which=fuse.which, K=K)
  S        <- Sgen(p=p, parg=parg, fuse.which=fuse.which, K=K, beta=beta.glm)
  B        <- Bgen(p=p, parg=parg, fuse.which=fuse.which, K=K, beta=S%*%beta.glm)
  
  theta.glm<- drop(B %*% S %*% beta.glm)   # cbind(beta.glm, solve(B%*%S)%*%theta.glm)
  X.t      <- X.b %*% solve(B%*%S)
  
  center   <- centerlist(p=p, parg=parg)
  al.r     <- 1
  al.weights <- 1/(abs(theta.glm)^al.r)
  al.weights[center] <- alpha * al.weights[center]  # if alpha!=0, there is sparsity penalty, except for intercept
  
  if(all(X.t[,1]==1)){
    glmnet.ada.fit <- glmnet(X.t[,-1], y, family=family, lambda=lambda, weights=weights, intercept=T, standardize=F, penalty.factor=al.weights[-1])
    thetahat       <- drop(coef(glmnet.ada.fit))
  } else{
    glmnet.ada.fit <- glmnet(X.t, y, family=family, lambda=lambda, weights=weights, intercept=F, standardize=F, penalty.factor=al.weights)
    thetahat       <- drop(coef(glmnet.ada.fit)[-1])
  }
  betahat   <- drop(solve(B%*%S)%*%thetahat)
  
  if.fuse       <- rep(0,p); if.fuse[fuse.which] <- 1; names(if.fuse) <- X.names
  out.beta      <- betahatReformat(beta=betahat, parg=parg, p=p, K=K); colnames(out.beta) <- X.names
  out.DF        <- as.matrix(betaDF(beta=betahat, parg=parg, p=p, X.names=X.names)); colnames(out.DF) <- "DF"
  
  out           <- list(family=family, alpha=alpha, if.fuse=if.fuse, betahat=out.beta, betainfo=t(out.DF))
  
  return(out)
}


X_trans <- function(X=X, fuse.which=fuse.which, K=K, X.names=X.names, sid=sid, studies=studies){
  p <- ncol(X)
  N <- nrow(X)
  X.b <- c()
  X.b.names <- c()
  for(i in 1:p){
    if(i %in% fuse.which){
      for(j in studies){
        temp <- rep(0,N)
        temp[(sid==j)] <- X[(sid==j),i]
        X.b <- cbind(X.b, temp)
      }
      X.b.names <- c(X.b.names, paste(X.names[i],"_",studies,sep=""))
    } else{
      X.b <- cbind(X.b, X[,i])
      X.b.names <- c(X.b.names, X.names[i])
    }
  }
  colnames(X.b) <- X.b.names
  return(X.b)
}

glm_wrap <- function(X=X, y=y, family=family, weights=weights){
  beta <- glm(y~.-1, data=data.frame(y,X), family=family, weights=weights)$coef
  return(beta)
}

pargroup <- function(p=p, fuse.which=fuse.which, K=K){
  parg <- c()
  for(i in 1:p){
    if(i %in% fuse.which){
      parg <- c(parg, rep(i,K))
    } else{
      parg <- c(parg, i)
    }
  }
  return(parg)
}

Bgen <- function(p=p, parg=parg, fuse.which=fuse.which, K=K, beta=beta){
  B <- c()
  for(i in 1:p){
    if(i %in% fuse.which){
      temp <- diag(c(0,rep(1,K-1)))
      temp[1,which.min(abs(beta[parg==i]))] <- 1
      for(i in 2:K){temp[i,i-1] <- -1}
      if(length(B)==0){B<-temp} else{B <- bdiag(B,temp)}
    } else{
      if(length(B)==0){B<-1} else{B <- bdiag(B,1)}
    }
  }
  return(as.matrix(B))
}

Sgen <- function(p=p, parg=parg, fuse.which=fuse.which, K=K, beta=beta){
  S <- c()
  for(i in 1:p){
    if(i %in% fuse.which){
      btemp <- beta[parg==i]
      order <- order(btemp)
      Stemp <- matrix(0,K,K)
      for(j in 1:K){
        Stemp[j,order[j]] <- 1
      }
      if(length(S)==0){S<-Stemp} else{S <- bdiag(S,Stemp)}
    } else{
      if(length(S)==0){S<-1} else{S <- bdiag(S,1)}
    }
  }
  return(as.matrix(S))
}

centerlist <- function(p=p, parg=parg){
  center <- sapply(c(1:p), function(x){ which(parg == x)[1] })
  return(center)
}

betaDF <- function(beta=beta, parg=parg, p=p, X.names=X.names){
  out <- c()
  for(i in 1:p){ out <- c(out, length(unique(beta[parg==i]))) }
  names(out) <- c(X.names)
  return(out)
}

critfunc <- function(betahat=betahat, family=family, criterion=criterion, N=N, X.b=X.b, y=y, K=K, p=p, sid=sid, parg=parg, fuse.which=fuse.which){
  mu <- 0; tau <- 0;
  for(i in c(1:p)){
    mu <- mu + length(unique(betahat[parg==i]))
    tau <- tau + choose(K, length(unique(betahat[parg==i])))
  }
  if(criterion == "AIC" ){ extra <- 2*mu }
  if(criterion == "BIC" ){ extra <- mu*log(N) }
  if(criterion == "EBIC"){ extra <- mu*log(N) + 2*1*log(tau)  }
  nk <- table(sid); nbar <- mean(nk)
  if(family=="gaussian") {LL <- sapply(c(1:K), function(kk) {-log(sum((y[sid==kk] - X.b[sid==kk,]%*%betahat)^2)/(nk[kk]))*nk[kk]/2})}
  if(family=="binomial"){LL <- sapply(c(1:K), function(kk) {t(y[sid==kk])%*%X.b[sid==kk,]%*%betahat - sum(log(1+exp(X.b[sid==kk,]%*%betahat)))})}
  if(family=="poisson") {LL <- sapply(c(1:K), function(kk) {t(y[sid==kk])%*%X.b[sid==kk,]%*%betahat - sum(exp(X.b[sid==kk,]%*%betahat))})}
  return(-2 * sum(LL * (nbar/nk)) + extra)
}

readkey <- function(){
  cat ("Press [enter] to continue")
  line <- readline()
}

plotCrit <- function(critPath=critPath, lamlist=lamlist, lambda_opt=lambda_opt, loglambda=F, criterion=criterion){
  xlabel <- expression(lambda)
  ylabel <- criterion
  if(loglambda==T) {lamlist <- log(lamlist+1, base=10)
                    lambda_opt <- log(lambda_opt+1, base=10)
                    xlabel <- expression(paste(log[10], "(", lambda, "+1)"))}
  plot(lamlist,critPath,type="l",pch=19,lwd=1.5, xlab=xlabel, ylab=ylabel, main="Model Selection")
  abline(v=lambda_opt, lty=3, lwd=1) # black dotted line
}

plotSolPath <- function(solPath=solPath, lamlist=lamlist, parg=parg, p=p, lambda_opt=lambda_opt, loglambda=F, X.names=X.names, fuse.which=fuse.which){
  sp <- solPath
  xlabel <- expression(lambda)
  if(loglambda==T) {lamlist <- log(lamlist+1, base=10)
                    lambda_opt <- log(lambda_opt+1, base=10)
                    xlabel <- expression(paste(log[10], "(", lambda, "+1)"))}
  plot(lamlist, sp[1,], type="n", ylim=range(sp), ylab=expression(beta), xlab=xlabel, main="Solution Path")
  allcolor <- rainbow(p+2)
  alltype  <- c(1:p)
  for(i in 1:p){
    # if(i==2) next
    for(j in which(parg==i)){
      # lines(lamlist,sp[j,],type="l", lwd=1.5, col=allcolor[i], lty=alltype[i]) # colored lines
      lines(lamlist,sp[j,],type="l", lwd=1, col=i) # colored lines
      # lines(lamlist,sp[j,], type="l", lwd=1.5, lty=alltype[i]) # dotted lines
    }
  }
  abline(v=lambda_opt, lty=3, lwd=1) # black dotted line
  # legend("topright", legend=X.names, lwd=1, col=allcolor, lty=alltype[i], cex=0.8)
  legend("topright", legend=X.names, lwd=1, col=c(1:p), cex=0.8)
  # legend("topright", legend=X.names, lwd=1, lty=alltype, cex=0.8)
}

betahatReformat <- function(beta=beta, parg=parg, p=p, K=K){
  beta_out <- list(); 
  for(i in 1:p){
    value <- beta[parg==i]
    if(length(value)==1){value <- rep(value,K)}
    group <- rep(NA, K)
    for(j in 1:length(unique(value))) group[value==unique(value)[j]] <- j
    listitem <- rbind(value, group)
    beta_out[[i]] <- listitem
  }
  beta_out_copy <- beta_out 
  beta_out <- c()
  for(i in 1:length(beta_out_copy)){beta_out <- cbind(beta_out, beta_out_copy[[i]][1,])}
  rownames(beta_out) <- paste("study", c(1:nrow(beta_out)), sep="")
  return(beta_out)
}

getGrouping <- function(betahat=betahat, parg=parg, p=p, fuse.which=fuse.which){
  out <- c()
  for(i in 1:p){
    beta <- betahat[parg==i]
    group <- rep(NA, length(beta))
    for(j in 1:length(unique(beta))) group[beta==unique(beta)[j]] <- j
    out <- c(out, group)
  }
  return(out)
}

fusiogram <- function(groupPath=groupPath, lambdaCP=lambdaCP, beta=beta, K=K, parg=parg,
                      plot.which=plot.which, lambda_opt=lambda_opt, loglambda=F, X.names=X.names){
  groupPath <- groupPath[,rev(1:ncol(groupPath))]
  if(loglambda==T) {lambdaCP <- log(lambdaCP+1, base=10)}
  lambdaCP <- rev(lambdaCP)
  temp <- rep(0, nrow(groupPath))
  distmat <- c()
  addvalue <- rep(max(lambdaCP), nrow(groupPath))
  for(i in 2:(ncol(groupPath)-1)){
    changeindex <- rep(F, nrow(groupPath))
    for(j in 1:length(unique(groupPath[,i]))){
      if(length(unique(groupPath[groupPath[,i]==unique(groupPath[,i])[j],i] - 
                         groupPath[groupPath[,i]==unique(groupPath[,i])[j],i+1])) != 1){
        changeindex <- (changeindex | groupPath[,i]==unique(groupPath[,i])[j])
      }
    }
    temp <- addvalue*changeindex
    distmat <- cbind(distmat, temp)
    addvalue[changeindex] <- lambdaCP[i]
    if(i == (ncol(groupPath)-1)){
      temp1 <- rep(0, nrow(groupPath))
      for(k in 1:length(unique(addvalue))){
        temp11 <- which(addvalue == unique(addvalue)[k])
        if(length(temp11) <= 2) {
          temp1[temp11[1]] <- addvalue[temp11[1]]
        } else {
          for(kk in 1:length(unique(groupPath[temp11,i]))){
            temp111 <- which(groupPath[,i] == unique(groupPath[temp11,i])[kk])
            temp1[temp111[1]] <- addvalue[temp111[1]]
          }
        }
      }    
      distmat <- cbind(distmat, temp1)
    }
  }
  rownames(distmat) <- names(beta)[parg==plot.which]
  ylabel <- expression(lambda)
  if(loglambda == T){
    lambda_opt <- log(lambda_opt+1, 10)
    ylabel <- expression(paste(log[10], "(", lambda, "+1)"))
  }
  plot(hclust(dist(distmat, method="maximum"), method="single"), hang=-1, main=paste(X.names[plot.which], sep=""),
       ylab=ylabel, xlab="", sub="", ylim=range(lambdaCP), lwd=1)
  abline(h=lambda_opt, lty=3, lwd=1)
}


# for evaluation, calculating entropy and purity
eval.cls <- function(pred, true){
  L <- length(unique(true))
  J <- length(unique(pred))
  tb <- table(pred,true)
  m <- sum(tb)
  mi <- rowSums(tb)
  mj <- colSums(tb)
  
  pij <- tb * (1/mi) %*% t(rep(1,L))
  ei <- apply(pij, 1, function(x) sum((x*log(x))[!is.na(x*log(x))]) )
  entropy <- -sum((mi/m) * ei)
  purity <- sum(apply(pij, 1, max) * (mi/m))
  return(list(entropy=entropy, purity=purity))
}