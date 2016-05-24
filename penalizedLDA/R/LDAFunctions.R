#library(flsa)

permute.rows <- function(x){
    dd <- dim(x)
      n <- dd[1]
      p <- dd[2]
      mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
      matrix(t(x)[order(mm)], n, p, byrow = T)
  }


balanced.folds <- function(y, nfolds = min(min(table(y)), 10)){
    totals <- table(y)
      fmax <- max(totals)
      nfolds <- min(nfolds, fmax)
      # makes no sense to have more folds than the max class size
      folds <- as.list(seq(nfolds))
      yids <- split(seq(y), y)
      # nice way to get the ids in a list, split by class
      ###Make a big matrix, with enough rows to get in all the folds per class
      bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
      for(i in seq(totals)) {
            bigmat[seq(totals[i]), i] <- sample(yids[[i]])
          }
      smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
      ### Now do a clever sort to mix up the NAs
      smallmat <- permute.rows(t(smallmat)) ### Now a clever unlisting
      x <- apply(smallmat, 2, function(x) x[!is.na(x)])
      if(is.matrix(x)){
            xlist <- list()
                for(i in 1:ncol(x)){
                        xlist[[i]] <- x[,i]
                      }
                return(xlist)
          }
      return(x)
  }






MakeYMat <- function(y){
  return(diag(length(unique(y)))[y,])
}

#MakeYMat <- function(y){
#  ymat <- matrix(0, nrow=length(y), ncol=length(unique(y)))
#  for(i in 1:ncol(ymat)) ymat[y==i,i] <- 1
#  return(ymat)
#}

MakeMeanVecs <- function(x,y){
  Y <- MakeYMat(y)
  return(solve(msqrt(t(Y)%*%Y))%*%t(Y)%*%x)
}

msqrt <- function(mat){
  if(sum((mat-t(mat))^2)>1e-8) stop("Msqrt function only works if mat is symmetric....")
  redo <- TRUE
  while(redo){
    eigenmat <- eigen(mat)
    d <- eigenmat$values
    d[abs(d)<(1e-12)] <- 0
    a <- eigenmat$vectors%*%diag(sqrt(d))%*%t(eigenmat$vectors)
    if(sum(is.na(a))==0) redo <- FALSE
    if(redo) print('did one loop') 
  }
  return(a)
}

soft <- function(mat,lam){
  return(sign(mat)*pmax(abs(mat)-lam, 0))
}


Penalty <- function(v,lambda,type,chrom, lambda2){
  if(type=="standard") return(lambda*sum(abs(v)))
  if(type=="ordered"){
    tots <- lambda*sum(abs(v))
    for(chr in sort(unique(chrom))){
      tots <- tots+lambda2*sum(abs(diff(v[chrom==chr])))
    }
    return(tots)
  }
}


PenalizedPCACrit <- function(x, P, v, lambda, d, type, chrom, lambda2){
  return(t(v)%*%t(x)%*%P%*%x%*%v-d*Penalty(v, lambda, type, chrom, lambda2))
}


PenalizedPCA <- function(x, lambda, K, type="standard",  chrom=NULL, lambda2=NULL, maxiter=30, trace=FALSE){
  # Notice that this K is the number of components desired, NOT the number of classes in the classification problem.
  # Here, x is (# of classes) \times p

  # The criterion is maximize_b (b' Sigmabet b) - P(b) s.t. b' Sigmawit b = 1
  # Where Sigmawit=I and where P(b) = lambda||b||_1 or P(b) = lambda||b||_1 + lambda2 ||b_i - b_{i-1}||_1
  # We take a MINORIZATION approach to this problem.
  if(type=="ordered" && is.null(chrom)) chrom <- rep(1, ncol(x))
  if(is.null(lambda2)) lambda2 <- lambda
  crits <-  NULL
  betas <- matrix(0, nrow=ncol(x), ncol=K)
  critslist <- list()
  for(k in 1:K){
    if(trace) cat("Starting on component ", k, fill=TRUE)
    if(k>1){
      svda <- svd(x%*%betas)
      u <- svda$u[,svda$d>(1e-10)]
      P <- diag(nrow(x)) - u%*%t(u)
    }
    if(k==1) P <- diag(nrow(x))
    svdx <- svd(t(x)%*%P)
    d <- svdx$d[1]^2
    beta <- svdx$u[,1]
    crits <- c(crits, PenalizedPCACrit(x, P, beta, lambda, d, type, chrom=chrom, lambda2))
    for(iter in 1:maxiter){
      if((length(crits)<4 || abs(crits[length(crits)]-crits[length(crits)-1])/max(1e-3, crits[length(crits)]) > (1e-6)) && sum(abs(beta))>0){
        if(trace) cat(iter,fill=FALSE)
        tmp <- (t(x)%*%P)%*%(x%*%beta)
        if(type=="standard") beta <- soft(tmp, d*lambda/2)
        if(type=="ordered"){
          for(chr in sort(unique(chrom))){
            beta[chrom==chr] <- as.numeric(flsa(tmp[chrom==chr],  d*lambda/2,  d*lambda2/2))
          }
        }
        beta <- beta/l2n(beta)
        beta[is.na(beta)] <- 0
        crits <- c(crits, PenalizedPCACrit(x, P, beta, lambda, d, type, chrom=chrom, lambda2))
      }
    }
    if(trace) cat(fill=TRUE)
    betas[,k] <- beta#cbind(betas, beta)
    critslist[[k]] <- crits
    if(min(diff(crits))<(-1e-6)) stop("min diff crits is too small!!!")
    crits <- NULL
  }
  return(list(v=betas, crits=as.vector(critslist)))
}
 




l2n <- function(vec){
  return(sqrt(sum(vec^2)))
}






diag.disc <-function(x, centroids, prior) {
  dd <- t(x) %*% centroids
  dd0 <- (rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
  scale(dd, as.numeric(dd0), FALSE) # this is -.5*||x_i - mu_k||^2+log(pi_k)
}




Classify <- function(xtr,xte,ytr,equalpriors=FALSE){ # I introduced unequal priors on 02/22/2010
  prior <- rep(1/length(unique(ytr)), length(unique(ytr)))
  if(!equalpriors){             
    for(k in 1:length(unique(ytr))) prior[k] <- mean(ytr==k)
  }
  # classify test obs to nearest training centroid.
  if(is.matrix(xtr) && ncol(xtr)>1){
    mus <- matrix(0, nrow=ncol(xtr), ncol=length(unique(ytr)))
    for(k in 1:length(unique(ytr))){
       mus[,k] <- apply(xtr[ytr==k,], 2, mean)
     }
  } else {
    mus <- matrix(NA, nrow=1, ncol=length(unique(ytr)))
    for(k in 1:length(unique(ytr))) mus[1,k] <- mean(xtr[ytr==k])
  }
  negdists <- diag.disc(t(xte), mus, prior)
  return(apply(negdists,1,which.max))
}  





wcsd <- function(vec, y){
  K <- length(unique(y))
  n <- length(vec)
  tots <- 0
  for(k in unique(y)){
    tots <- tots + sum((vec[y==k]-mean(vec[y==k]))^2)
  }
  return(sqrt(tots/n))
}


wcsd.matrix <- function(x,Y){
  n <- nrow(x)
  return(sqrt((1/n)*apply(((diag(n)-Y%*%diag(1/apply(Y,2,sum))%*%t(Y))%*%x)^2,2,sum)))
}


