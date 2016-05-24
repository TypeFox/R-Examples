
## Generates DNA profiles of n unrelated individuals for a locus
genTypeRec <- function(x,t,n,z=rep(0,lx <- length(x))){
  if(t==0){
    y <- sample(length(x),size=n*2,prob=x,replace=TRUE)
  }
  else{
    y <- sample(length(x),size=1,prob=x)
    z[y] <- z[y]+1
    while(length(y)<(n*2)){
      a <- sample(lx,size=1,prob=x*(1-t)+z*t)
      y <- c(y,a)
      z[a] <- z[a]+1
    }
  }
  data.frame(a1=y[nn <- seq(from=1,by=2,len=n)],a2=y[nn+1])
}

## Generates DNA profiles of n individuals. These are formed as n/2 pairs for
## relatives with a IDB-vector given by k. I.e. the profiles are mutually unrelated between pairs.
genRypeRec <- function(x,t,k,n,print=FALSE){
  if(n==0) return(NULL)
  y <- numeric(0)
  a <- c(0,0)
  if(is.list(x)){ z <- x[[2]]; x <- x[[1]]; lx <- length(x); }
  else z <- rep(0,lx <- length(x))
  if(t==0){
    while(length(y)<(n*2)){
      z <- rep(0,lx <- length(x))
      a <- sample(lx,size=2,prob=x,replace=TRUE)
      IBD <- sample(0:2,size=1,prob=k)
      if(IBD==2) b <- a
      else if(IBD==1) b <- c(a[sample(1:2,size=1)],sample(lx,size=1,prob=x))
      else if(IBD==0) b <- sample(lx,size=2,prob=x,replace=TRUE)
      y <- c(y,a,b)
    }
  }
  else{
    while(length(y)<(n*2)){
      a[1] <- sample(lx,size=1,prob=x*(1-t)+z*t)
      z[a[1]] <- z[a[1]]+1
      a[2] <- sample(lx,size=1,prob=x*(1-t)+z*t)
      z[a[2]] <- z[a[2]]+1
      y <- c(y,a)
      b <- gRec(a,x,z,t,k) ## relative's profile, allele probs, allele counts, theta, k-vector
      z[b[1]] <- z[b[1]]+1
      z[b[2]] <- z[b[2]]+1
      y <- c(y,b)
      if(print) message(paste("Pair ",format(length(y)/4,width=2),": R1=(",a[1],",",a[2],") R2=(",b[1],",",b[2],") x=",paste(z,sep=",",collapse=","),sep=""),appendLF=TRUE)
    }
  }
  data.frame(a1=y[nn <- seq(from=1,by=2,len=n)],a2=y[nn+1])
}

## Help function that creates profile of relative number 2 given the profile of relative 1.
gRec <- function(a,p,x,t,k){
  b <- c(0,0) ## profile of the relative
  IBD <- sample(x=0:2,size=1,prob=k)
  lx <- length(x)
  if(IBD==2) b <- a
  else if(IBD==1){
    b[1] <- a[sample(1:2,size=1,prob=c(0.5,0.5))]
##    x[b[1]] <- x[b[1]]+1
## CHANGED DUE TO THORE'S COMMENT: I.e. First relative (denoted a) is
## sampled "from posterior database". First allele of second relative
## (denoted b) is sampled from the same distribution and second allele
## from uniform distribution of a's alleles.
    b[2] <- sample(lx,size=1,prob=p*(1-t)+x*t)
  }
  else if(IBD==0){
    b[1] <- sample(lx,size=1,prob=p*(1-t)+x*t)
    x[b[1]] <- x[b[1]]+1
    b[2] <- sample(lx,size=1,prob=p*(1-t)+x*t)
  }
  b
}

## Wrapper function that combines the locus-specific databases
genUnrelatedRec <- function(probs,theta=0,n){
  db <- as.data.frame(do.call("cbind",lapply(probs,genTypeRec,t=theta,n=n)))
  names(db) <- paste(rep(names(probs),each=2),1:2,sep=".")
  cbind(id=1:n,db)
}

## Help function that computes the sufficient statistic for computing P(A_i | x^n) where x^n
## is the allele counts in a sample of n alleles - the sufficient statistic.
countAlleles <- function(x,p,start=1){
  x1 <- x[,y <- seq(from=start,to=ncol(x)-1,by=2)]
  x2 <- x[,y+1]
  names(x2) <- names(x1)
  lapply(rbind(x1,x2,unlist(lapply(p,length))),function(z) table(c(1:z[length(z)],z[-length(z)]))-1)
}

## A function for binding lists of equal lengths
lbind <- function(l1,l2){
  if(length(l1)!=length(l2)) stop("Lists are not of same length")
  ll <- lapply(l1,function(x) list(x,x))
  for(i in 1:length(ll)) ll[[i]][[2]] <- l2[[i]]
  ll
}

## A wrapper function that generates a database with relatives included.
genRelated <- function(probs,theta=0,n,rel=c("FS"=0,"C"=0,"PC"=0,"A"=0,"U"=0)){
  if(!all((rel%%2)==0)) stop("not all 'rel' counts are even numbers")
  names(rel) <- c("FS","C","PC","A","U")
  dbn <- paste(rep(names(probs),each=2),1:2,sep=".")
  ## Unrelated: (k0,k1,k2) = c(1,0,0)
  db <- as.data.frame(do.call("cbind",lapply(probs,genTypeRec,t=theta,n=rel[["U"]]))) 
  names(db) <- dbn
  ## Full-sibs: (k0,k1,k2) = c(1,2,1)/4
  z <- countAlleles(db,p=probs)
  if(nrow(db.fs <- as.data.frame(do.call("cbind",lapply(lbind(probs,z),genRypeRec,t=theta,k=c(1,2,1)/4,n=rel[["FS"]]))))>0){
    names(db.fs) <- dbn
    db <- rbind(db,db.fs)
  }
  ## First cousins: (k0,k1,k2) = c(3,1,0)/4
  z <- countAlleles(db,p=probs)
  if(nrow(db.c <- as.data.frame(do.call("cbind",lapply(lbind(probs,z),genRypeRec,t=theta,k=c(3,1,0)/4,n=rel[["C"]]))))>0){
    names(db.c) <- dbn
    db <- rbind(db,db.c)
  }
  ## Parent-child: (k0,k1,k2) = c(0,1,0)
  z <- countAlleles(db,p=probs)
  if(nrow(db.pc <- as.data.frame(do.call("cbind",lapply(lbind(probs,z),genRypeRec,t=theta,k=c(0,1,0),n=rel[["PC"]]))))>0){
    names(db.pc) <- dbn
    db <- rbind(db,db.pc)
  }
  ## Avuncular: (k0,k1,k2) = c(1,1,0)/2
  z <- countAlleles(db,p=probs)
  if(nrow(db.a <- as.data.frame(do.call("cbind",lapply(lbind(probs,z),genRypeRec,t=theta,k=c(1,1,0)/2,n=rel[["A"]]))))>0){
    names(db.a) <- dbn
    db <- rbind(db,db.a)
  }
  cbind(id=1:n,db)
}


## ## A wrapper function that generates a database with relatives included.
## genRelated <- function(probs,theta=0,n,rel=c("FS"=0,"C"=0,"PC"=0,"A"=0,"U"=0)){
##   if(!all((rel%%2)==0)) stop("not all 'rel' counts are even numbers")
##   names(rel) <- c("FS","C","PC","A","U")
##   dbn <- paste(rep(names(probs),each=2),1:2,sep=".")
##   odd <- seq(from=1,to=2*length(probs)-1,by=2)
##   if(rel[["U"]]<2) db.u <- NULL
##   else{
##     ## Unrelated: (k0,k1,k2) = c(1,0,0)
##     db.u <- as.data.frame(do.call("cbind",lapply(probs,genTypeRec,t=theta,n=rel[["U"]]))) 
##     names(db.u) <- dbn
##   }
##   db <- db.u
##   if(rel[["FS"]]<2) db.fs <- NULL
##   else {
##     ## Full-sibs: (k0,k1,k2) = c(1,2,1)/4
##     z <- countAlleles(db,p=probs)
##     db.fs <- as.data.frame(do.call("cbind",lapply(lbind(probs,z),genRypeRec,t=theta,k=c(1,2,1)/4,n=rel[["FS"]]))) 
##     names(db.fs) <- dbn
##   }
##   db <- rbind(db,db.fs)
##   if(rel[["C"]]<2) db.c <- NULL
##   else{
##     ## First cousins: (k0,k1,k2) = c(3,1,0)/4
##     z <- countAlleles(db,p=probs)
##     db.c <- as.data.frame(do.call("cbind",lapply(lbind(probs,z),genRypeRec,t=theta,k=c(3,1,0)/4,n=rel[["C"]]))) 
##     names(db.c) <- dbn
##   }
##   db <- rbind(db,db.c)
##   if(rel[["PC"]]<2) db.pc <- NULL
##   else{
##     ## Parent-child: (k0,k1,k2) = c(0,1,0)
##     z <- countAlleles(db,p=probs)
##     db.pc <- as.data.frame(do.call("cbind",lapply(lbind(probs,z),genRypeRec,t=theta,k=c(0,1,0),n=rel[["PC"]]))) 
##     names(db.pc) <- dbn
##   }
##   db <- rbind(db,db.pc)
##   if(rel[["A"]]<2) db.a <- NULL
##   else{
##     ## Avuncular: (k0,k1,k2) = c(1,1,0)/2
##     z <- countAlleles(db,p=probs)
##     db.a <- as.data.frame(do.call("cbind",lapply(lbind(probs,z),genRypeRec,t=theta,k=c(1,1,0)/2,n=rel[["A"]]))) 
##     names(db.a) <- dbn
##   }
##   db <- rbind(db,db.a)
##   cbind(id=1:n,db)
## }

dbSimulate <- function(probs,theta=0,n=1000,relatives=NULL){
  if(is.null(relatives)) return(genUnrelatedRec(probs,theta=theta,n=n))
  else{
    if(length(relatives)!=4)
      stop("'relatives' need to be a vector of length 4, where the elements gives the number of PAIRS of c(FULL-SIBLINGS, FIRST-COUSINS, PARENT-CHILD, AVUNCULAR)")
    if(n-sum(2*relatives)<0) stop("Too many pairs of relatives specified compared to the total number of profiles 'n'.")
    return(genRelated(probs,theta=theta,n=n,rel=c(relatives*2,n-sum(relatives*2))))
  }
}
