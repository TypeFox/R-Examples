# March 10 2009 - This function does sparse multiple CCA as described in Witten & Tibshirani (2009) extensions to sparse CCA method.


UpdateW <- function(xlist, i, K, sumabsthis, ws, type="standard", ws.final){
  tots <- 0
  for(j in (1:K)[-i]){
    diagmat <- (t(ws.final[[i]])%*%t(xlist[[i]]))%*%(xlist[[j]]%*%ws.final[[j]])
    diagmat[row(diagmat)!=col(diagmat)] <- 0
    tots <- tots + t(xlist[[i]])%*%(xlist[[j]]%*%ws[[j]]) - ws.final[[i]]%*%(diagmat%*%(t(ws.final[[j]])%*%ws[[j]]))
  }
  if(type=="standard"){
    sumabsthis <- BinarySearch(tots, sumabsthis)
    w <- soft(tots, sumabsthis)/l2n(soft(tots, sumabsthis))
  } else {
    tots <- as.numeric(tots)
    tots <- tots/mean(abs(tots)) 
    w <- FLSA(tots,lambda1=sumabsthis,lambda2=sumabsthis)[1,1,]
#    flsa.out <- diag.fused.lasso.new(tots,lam1=sumabsthis)
#    lam2ind <- which.min(abs(flsa.out$lam2-sumabsthis))
#    w <- flsa.out$coef[,lam2ind]
    w <- w/l2n(w)
    w[is.na(w)] <- 0
  }
  return(w)
}

GetCrit <- function(xlist, ws, K){
  crit <- 0
  for(i in 2:K){
    for(j in 1:(i-1)){
      crit <- crit + t(ws[[i]])%*%t(xlist[[i]])%*%xlist[[j]]%*%ws[[j]]
    }
  }
  return(crit)
}

GetCors <- function(xlist, ws, K){
  cors <- 0
  for(i in 2:K){
    for(j in 1:(i-1)){
      thiscor  <-  cor(xlist[[i]]%*%ws[[i]], xlist[[j]]%*%ws[[j]])
      if(is.na(thiscor)) thiscor <- 0
      cors <- cors + thiscor
    }
  }
  return(cors)
}


ftrans <- function(x){ return(.5*log((1+x)/(1-x))) }

MultiCCA.permute <- function(xlist, penalties=NULL, ws=NULL, type="standard", nperms=10, niter=3, trace=TRUE, standardize=TRUE){
  call <- match.call()
  K <- length(xlist)
  for(k in 1:K){
    if(ncol(xlist[[k]])<2) stop("Need at least 2 features in each data set!")
    if(standardize) xlist[[k]] <- scale(xlist[[k]], T, T)
  }
  if(length(type)==1) type <- rep(type, K) # If type is just a single element, expand to make a vector of length(xlist)
          # Or type can have standard/ordered for each elt of xlist
  if(length(type)!=K) stop("Type must be a vector of length 1, or length(xlist)")
  if(sum(type!="standard" & type!="ordered")>0) stop("Each element of type must be standard or ordered.")
  if(is.null(penalties)){
    if(sum(type=="ordered")==K) stop("Do not run MultiCCA.permute with only ordered data sets and penalties unspecified,
                                      since we only choose tuning the parameter via permutations when type='standard'.")
    penalties <- matrix(NA, nrow=K, ncol=10)
    for(k in 1:K){
      if(type[k]=="ordered"){
        lam <- ChooseLambda1Lambda2(svd(xlist[[k]])$v[,1])
        penalties[k,] <- lam
      } else {
        penalties[k,] <- pmax(seq(.1, .8, len=10)*sqrt(ncol(xlist[[k]])),1.1)
      }
    }
  }
  numnonzeros <- NULL
  if(!is.matrix(penalties)) penalties <- matrix(1,nrow=K,ncol=1)%*%matrix(penalties,nrow=1)
  permcors <- matrix(NA, nrow=nperms, ncol=ncol(penalties))
  cors <- numeric(ncol(penalties)) 
  for(i in 1:ncol(penalties)){
    out <- MultiCCA(xlist, penalty=penalties[,i], niter=niter, type=type, ws=ws, trace=trace)
    cors[i] <- GetCors(xlist, out$ws, K)
    numnonzeros <- c(numnonzeros, sum(out$numnonzeros))
    ws.init  <- out$ws.init
  }
  cat(fill=TRUE)
  for(j in 1:nperms){
    if(trace) cat("Permutation ", j, "of " , nperms ,fill=TRUE)
    xlistperm <- xlist
    for(k in 1:K){
      xlistperm[[k]] <- xlistperm[[k]][sample(1:nrow(xlistperm[[k]])),]
    }
    for(i in 1:ncol(penalties)){
      out <- MultiCCA(xlistperm, penalty=penalties[,i], niter=niter, type=type, ws=ws, trace=FALSE)
      permcors[j,i] <- GetCors(xlistperm, out$ws, K)
    }
  }
  pvals =zs =  NULL
  for(i in 1:ncol(penalties)){
    pvals <- c(pvals, mean(permcors[,i]>=cors[i]))
    zs <- c(zs, (cors[i]-mean(permcors[,i]))/(sd(permcors[,i])+.05))
  }
  if(trace) cat(fill=TRUE)
  out <- list(pvals=pvals, zstat=zs, bestpenalties=penalties[,which.max(zs)], cors=cors, corperms=permcors, numnonzeros=numnonzeros, ws.init=ws.init, call=call, penalties=penalties, type=type, nperms=nperms)
  class(out) <- "MultiCCA.permute"
  return(out)
}

MultiCCA <- function(xlist, penalty=NULL, ws=NULL, niter=25, type="standard", ncomponents=1, trace=TRUE, standardize=TRUE){
  for(i in 1:length(xlist)){
    if(ncol(xlist[[i]])<2) stop("Need at least 2 features in each data set.")
  }
  call <- match.call()
  K <- length(xlist)
  if(length(type)==1) type <- rep(type, K) # If type is just a single element, expand to make a vector of length(xlist)
          # Or type can have standard/ordered for each elt of xlist
  if(length(type)!=K) stop("Type must be a vector of length 1, or length(xlist)")
  if(sum(type!="standard" & type!="ordered")>0) stop("Each element of type must be standard or ordered.")
  for(k in 1:K){
    if(standardize) xlist[[k]] <- scale(xlist[[k]], T, T)
  }
  if(!is.null(ws)){
    makenull <- FALSE
    for(i in 1:K){
      if(ncol(ws[[i]])<ncomponents) makenull <- TRUE
    }
    if(makenull) ws <- NULL
  }
  if(is.null(ws)){
    ws <- list()
    for(i in 1:K) ws[[i]] <- matrix(svd(xlist[[i]])$v[,1:ncomponents], ncol=ncomponents)
  }
  if(is.null(penalty)){
    penalty <- rep(NA, K)
    penalty[type=="standard"] <- 4 # this is the default value of sumabs
    for(k in 1:K){
      if(type[k]=="ordered"){
        v <- svd(xlist[[k]])$v[,1]
        penalty[k] <- ChooseLambda1Lambda2(v)
      }
    }
  }
  ws.init <- ws
  if(length(penalty)==1) penalty <- rep(penalty, K)
  if(sum(penalty<1 & type=="standard")) stop("Cannot constrain sum of absolute values of weights to be less than 1.")
  for(i in 1:length(xlist)){
    if(type[i]=="standard" && penalty[i]>sqrt(ncol(xlist[[i]]))) stop("L1 bound of weights should be no more than sqrt of the number of columns of the corresponding data set.", fill=TRUE)
  }
  ws.final <- list()
  for(i in 1:length(ws)) ws.final[[i]] <- matrix(0, nrow=ncol(xlist[[i]]), ncol=ncomponents)
  cors <- NULL
  for(comp in 1:ncomponents){
    ws <- list()
    for(i in 1:length(ws.init)) ws[[i]] <- ws.init[[i]][,comp]
    curiter <- 1
    crit.old <- -10
    crit <- -20
    storecrits <- NULL
    while(curiter<=niter && abs(crit.old-crit)/abs(crit.old)>.001 && crit.old!=0){
      crit.old <- crit
      crit <- GetCrit(xlist, ws, K)
      storecrits <- c(storecrits,crit)
      if(trace) cat(curiter, fill=FALSE)
      curiter <- curiter+1
      for(i in 1:K){
        ws[[i]] <- UpdateW(xlist, i, K, penalty[i], ws, type[i], ws.final)
      }
    }
    for(i in 1:length(ws)) ws.final[[i]][,comp] <- ws[[i]]
    cors <- c(cors, GetCors(xlist, ws,K))
  }
  out <- list(ws=ws.final, ws.init=ws.init, K=K, call=call, type=type, penalty=penalty, cors=cors)
  class(out) <- "MultiCCA"
  return(out)
}

print.MultiCCA <- function(x,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Sum_{i<j} Cor(Xi wi, Xj wj) = ", sep="", paste(round(x$cors,4), sep="", " "),fill=TRUE)
  cat("There are ", x$K, " data sets.", fill=TRUE)
  for(i in 1:x$K){
    cat("Data set ", i, " is of type ", x$type[i], ".", fill=TRUE)
    if(x$type[i]=="ordered") cat("Tuning parameter used: Lambda was ", round(x$penalty[i],4),fill=TRUE)
    if(x$type[i]=="standard") cat("Tuning parameter used: Sum(abs(w)) was ", round(x$penalty[i],4),fill=TRUE)
    cat("Num non-zero elements of canonical variate(s) for data set ", i, ":    ")
    if(is.matrix(x$ws[[i]])) cat(apply(x$ws[[i]]!=0, 2, sum), fill=TRUE)
    if(!is.matrix(x$ws[[i]])) cat(sum(x$ws[[i]]!=0), fill=TRUE)
    cat(fill=TRUE)
  }
}




print.MultiCCA.permute <- function(x,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  tab <- round(cbind(x$pvals, x$zstat, x$cors, colMeans(x$corperms)), 3)
  dimnames(tab) <- list(paste("Tuning parameter set ", sep="", 1:length(x$pvals)), c("P-Value", "Z", "Cors", "Cors Perm"))
  print(tab, quote=FALSE)
  cat("Highest z score: ", max(x$zstat), "\n")
  cat("P-value corresponding to highest z score: ", x$pvals[which.max(x$zstat)], fill=TRUE)
  cat("Tuning parameters corresponding to highest z score: ", round(x$bestpenalties,3), "\n")
  sumabslamvecs <- round(x$penalties,4)
  dimnames(sumabslamvecs) <- list(paste(paste("Data set", sep=" ", 1:nrow(sumabslamvecs)),
                                        sep="", paste(paste("; Type is ", sep="", x$type),sep="",": ")), 1:ncol(sumabslamvecs))
  cat(fill=TRUE)
  cat("Tuning parameters used: ",fill=TRUE)
  print(sumabslamvecs,quote=FALSE)
}

plot.MultiCCA.permute <- function(x,...){
  sumabss <- x$penalties
  ccs <- x$cors
  nperms <- x$nperms
  zstats <- x$zstat
  ccperms <- x$corperms
  par(mfrow=c(2,1))
  plot(1:ncol(sumabss), ccs, main="Correlations For Real/Permuted Data", xlab="Index of Tuning Parameter Set",
       ylab="Correlations", ylim=range(ccperms,ccs))
  points(1:ncol(sumabss),ccs,type="l")
  for(i in 1:nperms){
    points(1:ncol(sumabss),ccperms[i,],col="green")
  }
  plot(1:ncol(sumabss),zstats,main="Z", xlab="Index of Tuning Parameter Set", ylab="Z score")
  lines(1:ncol(sumabss),zstats)
}
