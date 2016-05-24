# predict on new data using foba model, or return coefficients
#
"predict.foba" <- function(object, newx, k, type=c("fit","coefficients"), ...) {
  type <- match.arg(type)
  ii <- features.from.path(object$path,k)

  beta=object$beta[,ii$j]
  b=drop(object$meany-t(object$meanx)%*%beta)

  vn <- names(object$meanx)
  vs <- ii$path
  names(vs) <- vn[abs(vs)]
  
  if (type == "coefficients") {
    obj=list(coefficients=beta,intercept=b, selected.variables=vs)
    return (obj)
  }

  y= as.vector(as.matrix(newx)%*%beta+b,mode="numeric")

  obj=list(coefficients=beta,intercept=b, selected.variables=vs,fit=y)
  return (obj)
}

# print the variable selection path of a foba model
#
"print.foba" <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  vn <- names(x$meanx)
  cat(paste(x$type," variable selection path:\n"))

  path=rbind(matrix(x$path,nrow=1),matrix(1:length(x$path),nrow=1))
  dimnames(path) <- list(c("Var", "Step"), vn[abs(x$path)])
  print(path)
  invisible(x)
}

# forward/backward greedy steps
#
foba <-function(x,y, type=c("foba","foba.aggressive", "foba.conservative", "forward","backward"), steps=0, intercept=TRUE, nu=0.5,lambda=1e-10) {
  call <- match.call()

  type <- match.arg(type)
  TYPE <- switch(type,
                 foba = "FoBa",
                 foba.aggressive = "FoBa (aggressive)",
                 foba.conservative = "FoBa (conservative)",
                 forward = "Forward Greedy",
                 backward = "Backward Greedy")

  n <- dim(x)[1]
  p <- dim(x)[2]

  if (lambda<0) {
    lambda=0
  }
  
  # centering data
  if(intercept){
    meanx <- drop(colMeans(x))
    x <- scale(x, meanx, FALSE) # centers x
    meany <- mean(y)
    y <- drop(y - meany)
  }
  else {
    meanx <- rep(0,p)
    meany <- 0
    y <- drop(y)
  }
  xscale=sqrt(colSums(x*x)+lambda)
  x <- scale(x,FALSE,xscale)
  names(meanx) <- dimnames(x)[[2]]

  s=steps
  if (type=="forward") {
    #forward greedy algorithm
    if (s==0) s=p;

    path=rep(0,s);
    beta=matrix(rep(0,s*p),ncol=s)
    
    r=y;
    for (k in 1:s) {
      ik=which.max(abs(t(r)%*%x))
      if (length(which(path[1:k]==ik))>0) {
        path=path[1:(k-1)]
        beta=beta[,1:(k-1)]
        break;
      }

      path[k]=ik;
      myfs=as.matrix(x[,path[1:k]]);
      w=myridge(myfs,y,lambda);
      r=myfs%*%w-y;
      beta[path[1:k],k]=w/xscale[path[1:k]]
    }

    object <- list(call=call, type=TYPE, path=path, beta = beta, meanx = meanx, meany=meany)
    class(object) <- "foba"
    return (object)
  }

  if (type=="backward") {
    # backward greedy algorithm
    s=p;
    path=rep(0,s);
    beta=matrix(rep(0,s*p),ncol=s)

    stat=c(1:p)

    for (i in 1:p) {
      k=p+1-i
      ii=stat[1:k]
      myfs=as.matrix(x[,ii])
      w=myridge(myfs,y,lambda);
      beta[ii,k]=w/xscale[ii]

      ik=which.min(abs(w))
      path[k]=stat[ik];
      stat[ik]=stat[k];
    }

    object <- list(call=call, type=TYPE, path=path, beta = beta, meanx = meanx, meany=meany)
    class(object) <- "foba"
    return (object)
  }

  # now the FoBa algorithm
  #
  if (s==0) {
    s=2*p+1
  }
  
  path=rep(0,s)
  beta=matrix(rep(0,s*p),ncol=s)
  stat=matrix(rep(0,2*s),nrow=2);

  r=y
  v=(t(y)%*%y)
  k=0;

  if (nu<0) { nu=0}
  if (nu>0.99) {nu=0.99}
    
  it =0
  minfw=v
  while (it < s) {
    # forward step
    ik=which.max(abs(t(r)%*%x))
    if (length(which(stat[1,1:k]==ik))>0) {
      path=path[1:it]
      beta=beta[,1:it]
      break;
    }
    k=k+1;

    stat[1,k]=ik;
    myfs=as.matrix(x[,stat[1,1:k]]);

    w=myridge(myfs,y,lambda)

    rp=myfs%*%w-y;
    vp=(t(rp)%*%rp)+lambda*(t(w)%*%w)
    deltak=v-vp
    r=rp
    v=vp;
    
    stat[2,k]=deltak;

    it = it+1

    path[it]=ik;
    beta[stat[1,1:k],it]=w/xscale[stat[1,1:k]]

    if (minfw> deltak) {
      minfw=deltak
    }
    #backward step
    totfw=0;
    totbw=0;
    while ((k>1) & (it<s)) {
      ik=which.min(abs(w));

      tmp=stat[1,k];
      stat[1,k]=stat[1,ik];
      stat[1,ik]=tmp;

      myfs=as.matrix(x[,stat[1,1:(k-1)]]);
      w=myridge(myfs,y,lambda)
      rp=myfs%*%w-y;
      vp=(t(rp)%*%rp)+lambda*(t(w)%*%w)
      deltak=vp-v

      if (type == "foba") {
        # standard foba
        if (deltak >= stat[2,k]*nu) {
          break;
        }
      }
      else {
        if (type=="foba.aggressive") {
          # aggressive foba
          if ((totbw+deltak) >= (totfw+stat[2,k])*nu) {
            break;
          }
          totfw=totfw+stat[2,k]
          totbw=totbw+deltak
        }
        else {
           # conservative foba
          if (deltak >= minfw*nu) {
            break;
          }
        }
      }
      it=it+1
      
      path[it]=-stat[1,k];
      k=k-1
      beta[stat[1,1:k],it]=w/xscale[stat[1,1:k]]
      
      r=rp
      v=vp
    }
  }

  object <- list(call=call, type=TYPE, path=path, beta = beta, meanx = meanx, meany=meany)
  class(object) <- "foba"
  return (object)               

}

# a simple but inefficient ridge regression solver:
# a more efficient implementation requires keeping track of
# rank-one updates after each greedy step
#
myridge <- function(x,y,lambda) {
  x <- as.matrix(x)
  w=solve(t(x)%*%x + diag(lambda,dim(x)[2]),t(x)%*%y);
  return(w)
}
  
#
# best s features from FoBa path: last point in the path with s features
#
features.from.path <- function(path, s) {
  k=0;
  bestj=1
  besti=path[1]
  mys=1;
  kk=rep(0,s);
  for (j in 1:length(path)) {

    if (path[j]>0) {
      k=k+1;
      kk[k]=path[j];
    }
    else {
      ik=which(kk[1:k]==-path[j]);
      kk[ik[1]]=kk[k];
      k=k-1;
    }
    if ((mys<k) & (k<=s)) {
      mys=k
    }
    if (k==mys) {
      besti=kk[1:k];
      bestj=j;
    }
  }
  return(list(path=besti,j=bestj));
}


