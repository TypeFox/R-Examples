
#'@import combinat

#'@rdname CorrBin-internal

    mChoose <- function(n, rvec, log=FALSE){
      rlast <- n - sum(rvec)
      rveclong <- c(rvec, rlast)
      if (any(rveclong < 0)) return(0)
      
      res <- lgamma(n + 1) - sum(lgamma(rveclong + 1))
      if (log) res else exp(res)
    }

#' Estimate joint event probabilities for multinomial data
#'
#' An exchangeable multinomial distribution with \eqn{K+1} categories \eqn{O_1,\ldots,O_{K+1}}, can be
#' parameterized by the joint probabilities of events
#'\deqn{\tau_{r_1,\ldots,r_{K}|n} = P\big[X_1=\cdots=X_{r_1}=O_1,\ldots, X_{\sum_{i=1}^{K-1}r_i+1} =\cdots=X_{\sum_{i=1}^{K}r_i}=O_K\big] }{tau_{r_1,..,r_K|n} = P[X_1=...=X_{r_1}=O_1,..., X_{sum_{i=1}^{K-1}r_i+1} =...=X_{sum_{i=1}^{K}r_i}=O_K]}
#'where \eqn{r_i \geq 0} and \eqn{r_1+\cdots +r_K\leq n}.
#'The \code{jointprobs} function estimates these probabilities under various settings. 
#'Note that when some of the \eqn{r_i}'s equal zero, then no restriction on the number of outcomes of the 
#' corresponding type are imposed, so the resulting probabilities are marginal.
#'
#'@param cmdata a \code{CMData} object
#'@param type character string describing the desired type of estimate:
#' \itemize{
#'  \item{"averaged"}{ - averaged over the observed cluster-size distribution within each treatment}
#'  \item{"cluster"}{ - separately for each cluster size within each treatment}
#'  \item{"mc"}{ - assuming marginal compatibility, ie that \eqn{\tau} does not depend on the cluster-size}
#' }
#'@return a list with an array of estimates for each treatment. For a multinomial distribution with
#' \eqn{K+1} categories the arrays will have either \eqn{K+1} or {K} dimensions, depending on whether 
#' cluster-size specific estimates (\code{type="cluster"}) or pooled estimates 
#' (\code{type="averaged"} or \code{type="mc"}) are requested. For the cluster-size specific estimates #' the first dimension is the cluster-size. Each additional dimension is a possible outcome. 
#'
#'@seealso \code{\link{mc.est}} for estimating the distribution under marginal compatibility,
#'\code{\link{uniprobs}} and \code{\link{multi.corr}} for extracting the univariate marginal event
#'probabilities, and the within-multinomial correlations from the joint probabilities.
#'@examples
#'data(dehp)
#'# averaged over cluster-sizes
#'tau.ave <- jointprobs(dehp, type="ave")
#'# averaged P(X1=X2=O1, X3=O2) in the 1500 dose group
#'tau.ave[["1500"]]["2","1"]  # there are two type-1, and one type-2 outcome
#'
#'#plot P(X1=O1) - the marginal probability of a type-1 event over cluster-sizes
#'tau <- jointprobs(dehp, type="cluster")
#'ests <- as.data.frame(lapply(tau, function(x)x[,"1","0"]))
#'matplot(ests, type="b")
#'@export

jointprobs <- function(cmdata, type=c("averaged","cluster","mc")){
  type <- match.arg(type)
  
  
    nc <- attr(cmdata, "ncat")
    nrespvars <- paste("NResp", 1:nc, sep=".")
    M <- max(cmdata$ClusterSize)
  
  # multinomial lookup table
  mctab <- mChooseTable(M, nc, log=FALSE)
  
  res <- list()
  for (trt in levels(cmdata$Trt)){
    cm1 <- cmdata[cmdata$Trt==trt,]
    # observed freq lookup table
    atab <- array(0, dim=rep(M+1, nc))
    a.idx <- data.matrix(cm1[,nrespvars])
    atab[a.idx + 1] <- atab[a.idx + 1] + cm1$Freq
    
    if (type=="averaged"){
      Mn <- sum(cm1$Freq)
      
          
      res.trt <- array(NA, dim=rep(M+1, nc-1))
      dimnames(res.trt) <- rep.int(list(0:M), nc-1) 
      names(dimnames(res.trt)) <- paste("R", 1:(nc-1), sep="")
      # indices for possible values of r

         idx  <- hcube(rep( M +1,  nc-1 ))-1
          idxsum  <- rowSums(idx )
         idx  <- idx [ idxsum  <=  M , ,drop=FALSE]  #remove impossible indices
          idxsum  <-  idxsum [ idxsum  <=  M ]
      
      #indices for possible values of s 
      # (one more column than for r - ensures summation over all n's)

         sidx  <- hcube(rep( M +1,  nc ))-1
          sidxsum  <- rowSums(sidx )
         sidx  <- sidx [ sidxsum  <=  M , ,drop=FALSE]  #remove impossible indices
          sidxsum  <-  sidxsum [ sidxsum  <=  M ]
      
      for (i in 1:nrow(idx)){
        r <- idx[i,]
        s.idx <- which(sidxsum <= M-sum(r))
        lower.idx <- sidx[s.idx, , drop=FALSE]
        upper.idx <- lower.idx + rep(c(r,0), each=nrow(lower.idx))
        res.trt[rbind(r)+1] <- 
          sum(mctab[lower.idx+1] / mctab[upper.idx+1] * atab[upper.idx+1]) / Mn
      }
      
    } else if (type=="cluster") {
      Mn <- xtabs(Freq ~ factor(ClusterSize, levels=1:M), data=cm1) 
      
      res.trt <- array(NA, dim=c(M, rep(M+1, nc-1))) #first dimension is 'n'
      dimnames(res.trt) <- c(list(1:M), rep.int(list(0:M), nc-1)) 
      names(dimnames(res.trt)) <- c("N",paste("R", 1:(nc-1), sep=""))
      for (n in which(Mn > 0)){
        # indices for possible values of r
        
           idx  <- hcube(rep( n +1,  nc-1 ))-1
            idxsum  <- rowSums(idx )
           idx  <- idx [ idxsum  <=  n , ,drop=FALSE]  #remove impossible indices
            idxsum  <-  idxsum [ idxsum  <=  n ]
        
        for (i in 1:nrow(idx)){
          r <- idx[i,]
          s.idx <- which(idxsum <= n-sum(r))
          lower.idx <- idx[s.idx, , drop=FALSE]
          upper.idx <- lower.idx + rep(r, each=nrow(lower.idx))
          lower.idx <- cbind(lower.idx, n-sum(r)-idxsum[s.idx])   #add implied last column
          upper.idx <- cbind(upper.idx, n-sum(r)-idxsum[s.idx])   #add implied last column
          res.trt[cbind(n,rbind(r)+1)] <- 
            sum(mctab[lower.idx+1] / mctab[upper.idx+1] * atab[upper.idx+1]) / Mn[n]
        }
      }
      
    } else if (type=="mc") {
      
         pim <- mc.est.raw(cm1)[[1]]  #only one treatment group
         res.trt <- tau.from.pi(pim)
      
    }
    
    # append treatment-specific result to result list
    res.trt <- list(res.trt)
    names(res.trt) <- trt
    res <- c(res, res.trt) 
  }
  attr(res, "type") <- type
  res
}

#'@rdname mc.est
#'@method mc.est CMData
#'@export
#'@param eps numeric; EM iterations proceed until the sum of squared changes fall below \code{eps}  
#'@return For \code{CMData}: A data frame giving the estimated pdf for each treatment and
#'clustersize.  The probabilities add up to 1
#'for each \code{Trt}/\code{ClusterSize} combination. It has the following columns: 
#'@return \item{Prob}{numeric, the probability of \code{NResp} responses in a
#'cluster of size \code{ClusterSize} in group \code{Trt}}
#'@return \item{Trt}{factor, the treatment group}
#'@return \item{ClusterSize}{numeric, the cluster size}
#'@return \item{NResp.1 - NResp.K}{numeric, the number of responses of each type}
#'
#'@note
#'For multinomial data, the implementation is curerntly written in R, so it is not very fast.
#'
#'@examples
#'data(dehp)
#'dehp.mc <- mc.est(subset(dehp, Trt=="0"))
#'subset(dehp.mc, ClusterSize==2)

mc.est.CMData <- function(object, eps=1E-6, ...){

    nc <- attr(object, "ncat")      
    resp.vars1 <- paste("NResp", 1:(nc-1), sep=".")
   
    res <- mc.est.raw(object=object, eps=eps, ...)
    margres <- lapply(res, Marginals)  # has only NResp.1 - NResp.K
    
    mat.to.df <- function(idx, alist){
        dd <- as.data.frame.table(alist[[idx]], responseName="Prob")
        dd[c("N", resp.vars1)] <- lapply(dd[c("N", resp.vars1)], function(x)as.numeric(as.character(x)))
        dd$Trt <- names(alist)[idx]
        dd
    }
    margres <- lapply(1:length(margres), mat.to.df, alist=margres)
    fin <- do.call(rbind, margres)
    names(fin)[1] <- "ClusterSize"
    last.resp <- paste("NResp", nc, sep=".")
    fin[last.resp] <- fin$ClusterSize - rowSums(fin[resp.vars1]) # calculated omitted frequency
    fin$Trt <- factor(fin$Trt)
    fin <- fin[fin[last.resp] >= 0,]  #remove impossible clusters
    fin[c("Trt","ClusterSize", resp.vars1, last.resp, "Prob")]
}

#'@rdname CorrBin-internal
Marginals <- function(theta){
  K <- length(dim(theta))
  M <- dim(theta)[1]-1
  
  res <- array(0, dim=c(M, rep(M+1, K)))
  dimnames(res) <- c(N=list(1:M), dimnames(theta))
  
  # indices for possible values of r
  
     idx  <- hcube(rep( M +1,  K+1 ))-1
      clustersize  <- rowSums(idx )
     idx  <- idx [ clustersize  <=  M , ,drop=FALSE]  #remove impossible indices
      clustersize  <-  clustersize [ clustersize  <=  M ]
  
  idx <- idx[ , -1, drop=FALSE]  #remove (K+1)st category
  
  
    curridx <- idx[clustersize==M, ,drop=FALSE]
    res[cbind(M, curridx+1)] <- theta[curridx+1]
  
  for (cs in seq.int(M-1,1)){
    
      curridx <- idx[clustersize==cs, , drop=FALSE]
      res[cbind(cs, curridx+1)] <- (cs+1- rowSums(curridx))/(cs+1) * res[cbind(cs+1, curridx+1)]
      for (j in 1:K){
        lookidx <- curridx
        lookidx[ ,j] <- lookidx[ ,j] + 1   #add 1 to the j-th coordinate
        res[cbind(cs, curridx+1)] <- res[cbind(cs, curridx+1)] + 
                                     lookidx[,j]/(cs+1) * res[cbind(cs+1, lookidx+1)]
      }  
    
  }
  
  res
}

#'@rdname CorrBin-internal
mc.est.raw <- function(object, ...) UseMethod("mc.est.raw")

#'@method mc.est.raw CMData
mc.est.raw.CMData <- function(object, eps=1E-6, ...){
  cmdata <- object
  
    nc <- attr(cmdata, "ncat")
    nrespvars <- paste("NResp", 1:nc, sep=".")
    M <- max(cmdata$ClusterSize)
  
  
  # indices for possible values of r with clustersize = M
  
     idx  <- hcube(rep( M +1,  nc-1 ))-1
      idxsum  <- rowSums(idx )
     idx  <- idx [ idxsum  <=  M , ,drop=FALSE]  #remove impossible indices
      idxsum  <-  idxsum [ idxsum  <=  M ]
  

  res <- list()
  for (trt in levels(cmdata$Trt)){
    cm1 <- cmdata[cmdata$Trt==trt,]
    if (nrow(cm1) > 0){
      # observed freq lookup table
      atab <- array(0, dim=rep(M+1, nc))
      a.idx <- data.matrix(cm1[,nrespvars])
      atab[a.idx + 1] <- atab[a.idx + 1] + cm1$Freq
      Mn <- sum(cm1$Freq)
      
      
        res.trt <- array(NA, dim=rep(M+1, nc-1))
         
        #starting values
        res.trt[idx + 1] <- 1/nrow(idx)
        
        sqerror <- 1
        #EM update
        while (sqerror > eps){
              sqerror <- 0
              marg <- Marginals(res.trt)
          res.new <- array(NA, dim=rep(M+1, nc-1))
          res.new[idx + 1] <- 0
          
          
            for (i in 1:nrow(cm1)){
              rlong <- data.matrix(cm1[,nrespvars])[i,]    #nc elements
              r <- rlong[-nc]              #without the last category
              n <- cm1$ClusterSize[i]  
              # indices to which this cluster type contributes
              s.idx <- which(idxsum <= M-sum(r))
              tidx <- idx[s.idx, , drop=FALSE] + rep(r, each=length(s.idx))
              
              hvals <- apply(tidx, 1, function(tvec)prod(choose(tvec, r)) * choose(M-sum(tvec), n-sum(r))) 
              hvals <- hvals / choose(M, n)
              res.new[tidx+1] <- res.new[tidx+1] + atab[rbind(rlong)+1] / marg[rbind(c(n,r+1))] / Mn *
                                                   hvals * res.trt[tidx+1]
            }
          
              
          sqerror <- sum((res.new[idx+1] - res.trt[idx+1])^2)
              res.trt <- res.new 
        }
      
      
      # append treatment-specific result to result list
      dimnames(res.trt) <- rep.int(list(0:M), nc-1)
      names(dimnames(res.trt)) <- paste("NResp", 1:(nc-1), sep=".")
      res.trt <- list(res.trt)
    } else {
      res.trt <- list(c())
    } 
    res <- c(res, res.trt) 
  }
  names(res) <- levels(cmdata$Trt)
  res
} 
#'@rdname CorrBin-internal
tau.from.pi <- function(pimat){
  K <- length(dim(pimat))
  n <- dim(pimat)[1] - 1
  res <- array(NA, dim=rep(n+1, K)) 
  dimnames(res) <- rep.int(list(0:n), K) 
  names(dimnames(res)) <- paste("R", 1:K, sep="")

  # multinomial lookup table
  mctab <- mChooseTable(n, K+1, log=FALSE)
  
  # indices for possible values of r
  
     idx  <- hcube(rep( n +1,  K ))-1
      idxsum  <- rowSums(idx )
     idx  <- idx [ idxsum  <=  n , ,drop=FALSE]  #remove impossible indices
      idxsum  <-  idxsum [ idxsum  <=  n ]
  
  for (i in 1:nrow(idx)){
    r <- idx[i,]
    s.idx <- which(idxsum <= n-sum(r))
    lower.idx <- idx[s.idx, , drop=FALSE]
    upper.idx <- lower.idx + rep(r, each=nrow(lower.idx))
    lower.mc.idx <- cbind(lower.idx, n-sum(r)-idxsum[s.idx])   #add implied last column
    upper.mc.idx <- cbind(upper.idx, n-sum(r)-idxsum[s.idx])   #add implied last column
    res[rbind(r)+1] <- 
      sum(mctab[lower.mc.idx+1] / mctab[upper.mc.idx+1] * pimat[upper.idx+1])
  } 
  res
}

#'@rdname CorrBin-internal
p.from.tau <- function(taumat){
  K <- length(dim(taumat))
  idx <- diag(nrow=K)
  taumat[rbind(idx+1)]
}    

#'Extract univariate marginal probabilities from joint probability arrays
#'
#'Calculates the marginal probability of each event type for exchangeable correlated multinomial
#'data based on joint probability estimates calculated by the \code{\link{jointprobs}} function.
#'@param jp the output of \code{\link{jointprobs}} - a list of joint probability arrays by treatment
#'@param type one of c("averaged","cluster","mc") - the type of joint probability. By default,
#'the \code{type} attribute of \code{jp} is used.
#'@return a list of estimated probability of each outcome by treatment group. The elements are either
#'matrices or vectors depending on whether cluster-size specific estimates were requested 
#' (\code{type="cluster"}) or not.
#'@export
#'@seealso \code{\link{jointprobs}} for calculating the joint probability arrays
#'@examples
#'data(dehp)
#'tau <- jointprobs(dehp, type="averaged")
#'uniprobs(tau)
#'
#'#separately for each cluster size
#'tau2 <- jointprobs(dehp, type="cluster")
#'uniprobs(tau2)

uniprobs <- function(jp, type=attr(jp, "type")){
  type <- match.arg(type, c("averaged","cluster","mc"))
  
  get.probs <- function(tt){
    p <- p.from.tau(tt)
    c(p, 1-sum(p)) #add probability of last event type
  }

  if (type=="cluster") {
    res <- lapply(jp, function(x)apply(x, 1, get.probs))
  } else {
    res <- lapply(jp, get.probs)
  }
  res  
}

#'@rdname CorrBin-internal
corr.from.tau <- function(taumat){
  K <- length(dim(taumat))
  
  idx <- diag(nrow=K)
  numerator <- outer(1:K, 1:K, function(i,j){
     taumat[idx[i,]+idx[j,]+1] - taumat[idx[i,]+1] * taumat[idx[j,]+1]})
  denominator <- outer(1:K, 1:K, function(i,j){
     taumat[idx[i,]+1] * ifelse(i==j, 1-taumat[idx[i,]+1], -taumat[idx[j,]+1])})  
  res <- numerator / denominator    #the negative sign is in the denominator
  res
}


#'@rdname CorrBin-internal
corr.from.pi <- function(pimat){
  tt <- tau.from.pi(pimat)
  res <- corr.from.tau(tt)
  res
}


#'Extract correlation coefficients from joint probability arrays
#'
#'Calculates the within- and between-outcome correlation coefficients for exchangeable correlated
#'multinomial data based on joint probability estimates calculated by the \code{\link{jointprobs}}
#'function. These determine the variance inflation due the cluster structure.
#'
#'If \eqn{R_i} and \eqn{R_j} is the number of events of type \eqn{i} and \eqn{j}, respectively, in a cluster of
#'size \eqn{n}, then
#'\deqn{Var(R_i)= n p_i (1-p_i)(1 + (n-1)\phi_{ii})}
#'\deqn{Cov(R_i,R_j)= -n p_i p_j (1 + (n-1)\phi_{ij})}
#'where \eqn{p_i} and \eqn{p_j} are the marginal event probabilities and \eqn{\phi_{ij}} are the correlation
#' coefficients computed by \code{multi.corr}.
#'@param jp the output of \code{\link{jointprobs}} - a list of joint probability arrays by treatment
#'@param type one of c("averaged","cluster","mc") - the type of joint probability. By default,
#'the \code{type} attribute of \code{jp} is used.
#'@return a list of estimated correlation matrices by treatment group. If cluster-size specific 
#' estimates were requested (\code{(type="cluster")}), then each list elements are a list of
#' these matrices for each cluster size.
#'@export
#'@seealso \code{\link{jointprobs}} for calculating the joint probability arrays
#'@examples
#'data(dehp)
#'tau <- jointprobs(dehp, type="averaged")
#'multi.corr(tau)
#'

multi.corr <- function(jp, type=attr(jp, "type")){
  type <- match.arg(type, c("averaged","cluster","mc"))
  
  if (type=="cluster") {
    K <- length(dim(jp[[1]])) - 1
    resmat <- lapply(jp, function(x)apply(x, 1, corr.from.tau))
    res <- lapply(resmat, function(x){
                  lapply(1:ncol(x), function(idx)matrix(x[,idx], nrow=K))})
  } else {
    res <- lapply(jp, corr.from.tau)
  }
  res  
}

#'@rdname mc.test.chisq
#'@method mc.test.chisq CMData
#'@export
#'@examples
#'
#'data(dehp)
#'mc.test.chisq(dehp)
#' 

mc.test.chisq.CMData <- function(object, ...){
  cmdata <- object[object$Freq > 0, ]
  K <- attr(object, "ncat")-1
  nrespvars <- paste("NResp", 1:K, sep=".")
  
  get.T <- function(x){
      x$Trt <- factor(x$Trt)  #remove unused levels
      tt <- jointprobs(x, type="mc")[[1]] #only one treatment group
      p <- p.from.tau(tt)
      phi <- corr.from.tau(tt)
      xx <- x[rep(1:nrow(x), x$Freq),]
      xx$Freq <- 1
      
      M <- max(x$ClusterSize)
      Mn <- table(factor(xx$ClusterSize, levels=1:M)) 

      scores <- (1:M) - (M+1)/2
      
      Rmat <- data.matrix(xx[,nrespvars,drop=FALSE])
      nvec <- xx$ClusterSize
      cvec <- scores[nvec] 
      c.bar <- weighted.mean(cvec, w=nvec)
      cvec <- cvec - c.bar 
            
      X <- t(Rmat) %*% cvec
      Sigma <- diag(p, nrow=length(p)) - outer(p,p)  #multinomial vcov
      od.matrix <- matrix(0, nrow=K, ncol=K)  #over-dispersion matrix
      for (n in 1:M){
        od.matrix <- od.matrix + n * Mn[n] * (scores[n]-c.bar)^2 * (1+(n-1)*phi)
      }
      Sigma <- Sigma * od.matrix
      
      Tstat <- t(X) %*% solve(Sigma) %*% X       
      Tstat
   }
      
   chis <- by(cmdata, cmdata$Trt, get.T)
   chis <- chis[1:length(chis)]
   chi.list <- list(chi.sq=chis, p=pchisq(chis, df=K, lower.tail=FALSE))
   overall.chi <- sum(chis)
   overall.df <- length(chis) * K
   list(overall.chi=overall.chi, overall.p=pchisq(overall.chi, df=overall.df, lower.tail=FALSE), 
        individual=chi.list)
}

#'@rdname CorrBin-internal
  mChooseTable <- function(n, k, log=FALSE){
    res <- array(NA, dim=rep.int(n+1, k))
    dimnames(res) <- rep.int(list(0:n), k)
    
    idx <- hcube(rep.int(n+1, k)) - 1
    idx <- idx[rowSums(idx) <= n, ,drop=FALSE]
    for (i in 1:nrow(idx)){
        r <- idx[i, ]
        res[rbind(r)+1] <- mChoose(n=sum(r), rvec=r, log=log)
    }
    res
  }
