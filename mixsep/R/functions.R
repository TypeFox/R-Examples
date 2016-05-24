#### Functions used by the main functions

convertTab <- function(x,mode="paste",subset=NULL){
  if(!is.null(subset)){
    if(is.numeric(subset)) I <- (1:ncol(x))[subset]
    else if(is.character(subset)){
      I <- numeric(0)
      for(n in subset) I <- c(I,grep(n,names(x)))
    }
    else stop("Wrong use of subsetting")
  }
  else I <- 1:ncol(x)
  if(mode==substr("paste",0,nchar(mode))){ for(i in I) x[,i] <- paste(x[,i]) }
  if(mode==substr("numeric",0,nchar(mode))){ for(i in I) x[,i] <- as.numeric(x[,i]) }
  x
}

listsum <- function(x){
  if(length(x)==1) x[[1]]
  else x[[1]]+listsum(x[-1])
}

## Functions for constructing the covariance matrix
Cs <- function(n) diag(n)-matrix(1/n,n,n)

CA <- function(h){
  if(length(h)==1) Cs(length(h))%*%as.matrix(h)%*%t(Cs(length(h)))
  else Cs(length(h))%*%diag(h)%*%t(Cs(length(h)))
}

## MLE functions for alpha and tau under the covariance model
ahat.locus <- function(x){
  if(any(is.na(x$P1))) return(list(num=0,den=0))
  x0 <- (x$P1-x$P2)*sum(x$area)/2
  x1 <- x$P2*sum(x$area)/2
  iW <- ginv(CA(x$height))
  list(num=t(x0)%*%iW%*%(x$area-x1),den=t(x0)%*%iW%*%x0)
}

ahat <- function(x){
  numden <- lapply(x,ahat.locus)
  sum(unlist(lapply(numden,function(x) x$num)))/sum(unlist(lapply(numden,function(x) x$den)))
}

that.locus <- function(x,alpha){
  if(any(is.na(x$P1))) return(0)  
  if(missing(alpha)){ alpha <- ahat.locus(x); alpha <- alpha$num/alpha$den }
  x0 <- (x$P1-x$P2)*sum(x$area)/2
  x1 <- x$P2*sum(x$area)/2
  iW <- ginv(CA(x$height))
  t(x$area-alpha*x0-x1)%*%iW%*%(x$area-alpha*x0-x1)
}

that <- function(x,alpha){
  if(missing(alpha)) alpha <- ahat(x)
  Ds <- unlist(lapply(x,that.locus,alpha=alpha))
  sum(Ds)/(sum(unlist(lapply(x,function(y) (!any(is.na(y$P1)))*(nrow(y)-1))))-1)
}

## MLE functions for alpha and tau under the precision model 
atilde.locus <- function(x){
  if(any(is.na(x$P1))) return(list(num=0,den=0))
  x0 <- (x$P1-x$P2)*sum(x$area)/2
  x1 <- x$P2*sum(x$area)/2
  list(num=sum(x0*(x$area-x1)/x$height),den=sum(x0^2/x$height))
}

atilde <- function(x){
  numden <- lapply(x,atilde.locus)
  sum(unlist(lapply(numden,function(x) x$num)))/sum(unlist(lapply(numden,function(x) x$den)))
}

ttilde.locus <- function(x,alpha){
  if(missing(alpha)){ alpha <- atilde.locus(x); alpha <- alpha$num/alpha$den }
  if(any(is.na(x$P1))) return(0)
  x0 <- (x$P1-x$P2)*sum(x$area)/2
  x1 <- x$P2*sum(x$area)/2
  sum((x$area-alpha*x0-x1)^2/x$height)
}

ttilde <- function(x,alpha){
  if(missing(alpha)) alpha <- atilde(x)
  Ds <- unlist(lapply(x,ttilde.locus,alpha=alpha))
  sum(Ds)/(sum(unlist(lapply(x,function(y) (!any(is.na(y$P1)))*(nrow(y)-1))))-1)
}

## Function for constructing the profile list returned by mixsep()
mmm <- function(x){
  if(!is.element("P3",names(x))){ ## two-person mixture
    mm <- paste(sort(rep(paste(x$allele)[x$P1!=0],sum(x$P1))[1:2]),collapse=",")
    MM <- paste(sort(rep(paste(x$allele)[x$P2!=0],sum(x$P2))[1:2]),collapse=",")
    res <- matrix(c(mm,MM),2,1,dimnames=list("Profiles"=c("Minor","Major"),"Locus"=paste(x$sys)[1]))
  }
  if(is.element("P3",names(x))){ ## three-person mixture
    mm <- paste(sort(rep(paste(x$allele)[x$P1!=0],sum(x$P1))[1:2]),collapse=",")
    Mm <- paste(sort(rep(paste(x$allele)[x$P2!=0],sum(x$P2))[1:2]),collapse=",")
    MM <- paste(sort(rep(paste(x$allele)[x$P3!=0],sum(x$P3))[1:2]),collapse=",")
    res <- matrix(c(mm,Mm,MM),3,1,dimnames=list("Profiles"=c("Minor","Mid","Major"),"Locus"=paste(x$sys)[1]))
  }
  res
} 



############ m may be 2 or 3:

ahat.locus <- function(x,m=2){
  if(m==2){
    if(any(is.na(x$P1))) return(list(num=0,den=0))
    x0 <- (x$P1-x$P2)*sum(x$area)/2
    x1 <- x$P2*sum(x$area)/2
    iW <- ginv(CA(x$height))
    res <- list(num=t(x0)%*%iW%*%(x$area-x1),den=t(x0)%*%iW%*%x0)
  }
  if(m==3){
    if(any(is.na(x$P1))) return(list(num=c(0,0),den=c(0,0)))
    X <- matrix(c(x$P1-x$P3,x$P2-x$P3),nrow(x),2)*sum(x$area)/2
    Xm <- x$P3*sum(x$area)/2
    iW <- ginv(CA(x$height))
    res <- list(num=t(X)%*%iW%*%(x$area-Xm),den=t(X)%*%iW%*%X)
  }
  res
}

ahat <- function(x,m=2){
  numden <- lapply(x,ahat.locus,m=m)
  if(m==2) res <- sum(unlist(lapply(numden,function(x) x$num)))/sum(unlist(lapply(numden,function(x) x$den)))
  if(m==3) res <- ginv(listsum(lapply(numden,function(x) x$den)))%*%listsum(lapply(numden,function(x) x$num))
  res
}

that.locus <- function(x,alpha,m=2){
  if(m==2){
    if(any(is.na(x$P1))) return(0) 
    if(missing(alpha)){ alpha <- ahat.locus(x,m=m); alpha <- alpha$num/alpha$den }
    alpha <- ifelse(is.na(alpha),0,alpha)
    x0 <- (x$P1-x$P2)*sum(x$area)/2
    x1 <- x$P2*sum(x$area)/2
    iW <- ginv(CA(x$height))
    res <- t(x$area-alpha*x0-x1)%*%iW%*%(x$area-alpha*x0-x1)
  }
  if(m==3){
    if(any(is.na(x$P1))) return(c(0,0)) 
    if(missing(alpha)){ alpha <- ahat.locus(x,m=m); alpha <- solve(alpha$den)%*%alpha$num }
    alpha <- ifelse(is.na(alpha),0,alpha)
    X <- matrix(c(x$P1-x$P3,x$P2-x$P3),nrow(x),2)*sum(x$area)/2
    Xm <- x$P3*sum(x$area)/2
    iW <- ginv(CA(x$height))
    res <- t(x$area-X%*%alpha-Xm)%*%iW%*%(x$area-X%*%alpha-Xm)
  }
  res
}

that <- function(x,alpha,m=m){
  if(missing(alpha)) alpha <- ahat(x,m=m)
  Ds <- unlist(lapply(x,that.locus,alpha=alpha,m=m))
  sum(Ds)/(sum(unlist(lapply(x,function(y) (!any(is.na(y$P1)))*(nrow(y)-1))))-1)
}

## Tilde

atilde.locus <- function(x,m=2){
  if(m==2){
    if(any(is.na(x$P1))) return(list(num=0,den=0))
    x0 <- (x$P1-x$P2)*sum(x$area)/2
    x1 <- x$P2*sum(x$area)/2
    res <- list(num=sum(x0*(x$area-x1)/x$height),den=sum(x0^2/x$height))
  }
  if(m==3){
    if(any(is.na(x$P1))) return(list(num=c(0,0),den=c(0,0)))
    X <- matrix(c(x$P1-x$P3,x$P2-x$P3),nrow(x),2)*sum(x$area)/2
    Xm <- x$P3*sum(x$area)/2
    res <- list(num=t(X)%*%CA(1/x$height)%*%(x$area-Xm),den=t(X)%*%CA(1/x$height)%*%X)
  }
  res
}

atilde <- function(x,m=m){
  numden <- lapply(x,atilde.locus,m=m)
  if(m==2) res <- sum(unlist(lapply(numden,function(x) x$num)))/sum(unlist(lapply(numden,function(x) x$den)))
  if(m==3) res <- solve(listsum(lapply(numden,function(x) x$den)))%*%listsum(lapply(numden,function(x) x$num))
  res
}

ttilde.locus <- function(x,alpha,m=2){
  if(m==2){
    if(missing(alpha)){ alpha <- atilde.locus(x,m=m); alpha <- alpha$num/alpha$den }
    if(any(is.na(x$P1))) return(0)
    x0 <- (x$P1-x$P2)*sum(x$area)/2
    x1 <- x$P2*sum(x$area)/2
    res <- sum((x$area-alpha*x0-x1)^2/x$height)
  }
  if(m==3){
    if(missing(alpha)){ alpha <- atilde.locus(x,m=m); alpha <- solve(alpha$den)%*%alpha$num }
    if(any(is.na(x$P1))) return(c(0,0))
    X <- matrix(c(x$P1-x$P3,x$P3-x$P3),nrow(x),2)*sum(x$area)/2
    Xm <- x$P3*sum(x$area)/2
    res <- t(x$area-X%*%alpha-Xm)%*%CA(1/x$height)%*%(x$area-X%*%alpha-Xm)
  }
  res
}

ttilde <- function(x,alpha,m=m){
  if(missing(alpha)) alpha <- atilde(x,m=m)
  Ds <- unlist(lapply(x,ttilde.locus,alpha=alpha,m=m))
  sum(Ds)/(sum(unlist(lapply(x,function(y) (!any(is.na(y$P1)))*(nrow(y)-1))))-1)
}

locus.expectedAreas <- function(x,alpha,m=2){
  PP <- paste("P",1:m,sep="")
  x$exp <- rowSums(x[,PP]*matrix(rep(c(alpha,1-sum(alpha)),each=nrow(x)),nrow(x),m))*sum(x$area)/2
  x
}
    
expectedAreas <- function(x,alpha,m=2){
  do.call("rbind",lapply(x,locus.expectedAreas,alpha=alpha,m=m))
}

  
  
