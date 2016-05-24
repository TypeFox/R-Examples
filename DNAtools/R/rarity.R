
## Collapses a (m,p)-matrix to a (2*m+p)-vector
mpcollapse <- function(mallele,nloci){
  tmp <- data.frame(m2 = (m <- floor(mallele/2):0),m1=mallele-2*m)
  tmp[tmp$m1<=nloci & tmp$m2<=nloci,]
}

dbCollapse <- function(x){
  if(class(x)=="dbcompare") mpmatrix <- x$m
  else if(class(x)!="matrix") stop("Input must be of class 'matrix' or 'dbcompare'")
  else mpmatrix <- x
  res <- sapply(1:((nL <- (ncol(mpmatrix)-1))*2+1),function(i) sum(diag(mpmatrix[(MP <- mpcollapse(i-1,nL))$m2+1,,drop=FALSE][,MP$m1+1,drop=FALSE]),na.rm=TRUE))
  names(res) <- 0:(length(res)-1)
  res
}

wildCollapse <- function(x){
  maxAlleles <- nrow(x)-1
  wildCol <- rep(NA,maxAlleles+1)
  for(s in 0:maxAlleles) wildCol[s+1] <- sum(diag(as.matrix(x[0:s+1,rev(0:s+1)])))
  names(wildCol) <- 0:maxAlleles
  wildCol
}

## Extracts the upper left triangle of a quadratic matrix
up.tri <- function(x,diag=TRUE,droplast=FALSE){
  x <- as.matrix(x)
  res <- (row(x) + col(x)) - 1 <= ncol(x)
  if(!diag) res[(row(x) + col(x)-1) == ncol(x)] <- FALSE
  if(droplast) res[nrow(res),1] <- FALSE
  res
}

## Creates a matrix/list of cell names
dbCats <- function(nloci,vector=FALSE){
  res <- outer(0:nloci,0:nloci,function(i,j) paste(i,j,sep="/"))
  res[!up.tri(res)] <- NA
  if(vector) res <- t(res)[up.tri(res)]
  res
}

## The recursive step of the expected value function. See eq. (2) in Tvedebrink et al.
rare <- function(q){
  S <- ncol(q)
  M <- replicate(S,matrix(0,S+1,S+1),simplify=FALSE)
  M[[1]][1,1] <- q[1,1] # P_{0,1}
  M[[1]][1,2] <- q[2,1] # P_{1,1}
  M[[1]][2,1] <- q[3,1] # P_{2,1}
  for(s in 2:S){
    for(m in 1:(s+1)){
      for(p in 1:(s-m+2)){
        if(m==1 & p>1) M[[s]][m,p] <- q[1,s]*M[[s-1]][m,p]+q[2,s]*M[[s-1]][m,p-1]
        else if(m>1 & p==1) M[[s]][m,p] <- q[1,s]*M[[s-1]][m,p]+q[3,s]*M[[s-1]][m-1,p]
        else if(m==1 & p==1) M[[s]][m,p] <- q[1,s]*M[[s-1]][m,p]
        else M[[s]][m,p] <- q[1,s]*M[[s-1]][m,p]+q[2,s]*M[[s-1]][m,p-1]+q[3,s]*M[[s-1]][m-1,p]
      }
    }
  }
  x <- M[[S]]
  dimnames(x) <- list(0:S,0:S)
  x
}

## Computes the P_{0/0}, P_{0/1}, P_{1/0} for a given locus
Ps <- function(p,t,k=rep(0,3),r=0,R=0){ ## R and r not used
  s1 <- sum(p); s2 <- sum(p^2); s3 <- sum(p^3); s4 <- sum(p^4)
  d <- (1+t)*(1+2*t)
  p0 <- t^2*(1-t)*(1-s2) + 2*t*(1-t)^2*(1-2*s2+s3) + (1-t)^3*(1-4*s2+4*s3+2*s2^2-3*s4)
  p1 <- 8*t^2*(1-t)*(1-s2) + 4*t*(1-t)^2*(1-s3) + 4*(1-t)^3*(s2-s3-s2^2+s4)
  p2 <- 6*t^3 + t^2*(1-t)*(2+9*s2) + 2*t*(1-t)^2*(2*s2+s3) + (1-t)^3*(2*s2^2-s4)
  if(all(k==0)) res <- c(p0,p1,p2)/d
  else res <- c(k[3]*p0/d,k[2]*(1-t)*(1-s2)+k[3]*p1/d,k[1]+k[2]*(t+(1-t)*s2)+k[3]*p2/d)
  res
}

## The recursive step of the expected value function. See eq. (2) in Tvedebrink et al.
Frare <- function(q){
  S <- ncol(q)
  allCats <- as.vector(t(outer(0:S,0:S,FUN=paste,sep="/")))
  M <- replicate(S,matrix(0,(S+1)^2,(S+1)^2,dimnames=list(Genuine=paste(1:((S+1)^2),allCats,sep=":"),Wildcard=paste(1:((S+1)^2),allCats,sep=":"))),simplify=FALSE)
  M[[1]][1*(S+1)+(1+0),0*(S+1)+(1+0)] <- q[1,1] # P_{1/0/0/0} {m/p/Fm/Fp}
  M[[1]][0*(S+1)+(1+1),0*(S+1)+(1+0)] <- q[2,1] # P_{0/1/0/0}
  M[[1]][0*(S+1)+(1+0),1*(S+1)+(1+0)] <- q[3,1] # P_{0/0/1/0}
  M[[1]][0*(S+1)+(1+0),0*(S+1)+(1+1)] <- q[4,1] # P_{0/0/0/1}
  M[[1]][0*(S+1)+(1+0),0*(S+1)+(1+0)] <- q[5,1] # P_{0/0/0/0}
  for(s in 2:S){
    for(m in 1:(s+1)){
      for(p in 1:(s-m+2)){
        for(fm in 1:(s+1)){
          for(fp in 1:(s-fm+2)){
            if(m>1)  M[[s]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp] <- M[[s]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp] + q[1,s]*M[[s-1]][(m-2)*(S+1)+p,(fm-1)*(S+1)+fp]
            if(p>1)  M[[s]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp] <- M[[s]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp] + q[2,s]*M[[s-1]][(m-1)*(S+1)+(p-1),(fm-1)*(S+1)+fp]
            if(fm>1) M[[s]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp] <- M[[s]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp] + q[3,s]*M[[s-1]][(m-1)*(S+1)+p,(fm-2)*(S+1)+fp]
            if(fp>1) M[[s]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp] <- M[[s]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp] + q[4,s]*M[[s-1]][(m-1)*(S+1)+p,(fm-1)*(S+1)+(fp-1)]
            M[[s]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp] <- M[[s]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp] + q[5,s]*M[[s-1]][(m-1)*(S+1)+p,(fm-1)*(S+1)+fp]
          }
        }
      }
    }
  }
  x <- M[[S]]
  remCats <- as.vector(t(outer(0:S,0:S,FUN=function(x,y,n=S) (x+y)<=n)))
  dimnames(x) <- list(Genuine=allCats,Wildcard=allCats)
  x[remCats,remCats]
}

## Computes the P_{m/p/Fm/Fp}  for a given locus
FPs <- function(p,t,k=rep(0,3),r=0,R=0){ ## R and r not used
  s2 <- sum(p^2)
  s3 <- sum(p^3)
  s4 <- sum(p^4)
  d <- (1+2*t)*(1+t)
  pmatch <- 2*(t^2*(1-t)*(1-s2)+2*t*(1-t)^2*(s2-s3)+(1-t)^3*(s2^2-s4))
  pfmatch <- 6*t^3+t^2*(1-t)*(8+3*s2)+t*(1-t)^2*(12*s2-6*s3)+(1-t)^3*(4*s3-3*s4)
  p2 <- pmatch+pfmatch
  ppartial <- 4*(t*(1-t)^2*(1-3*s2+2*s3)+(1-t)^3*(s2-2*s3+2*s4-s2^2))
  pfpartial <- t^2*(1-t)*(1-s2)+t*(1-t)^2*(2-4*s2+2*s3)+(1-t)^3*(2*s2-4*s3+3*s4-1*s2^2)
  p1 <- ppartial+pfpartial
  pmismatch <- (1-t)^3*(1-6*s2+8*s3-6*s4+3*s2^2)
  p0 <- pmismatch
  if(all(k==0)) res <- c(p0,p1,p2)/d
  else res <- c(k[3]*p0/d,k[2]*(1-t)^2*(1-3*s2+2*s3)/(1+t)+k[3]*p1/d,
                k[1]+k[2]*(2*t^2+3*t*(1-t)+(1-t)^2*(3*s2-2*s3))/(1+t)+k[3]*p2/d)
  res
}

## Computes the P_{m/p/Fm/Fp}  for a given locus
RPs <- function(p,t,r=0,R=0,k=rep(0,3)){
    pp <- prob("all",p,R,r,t) ## Calls c++ functions for computing the relevant probabilities
    ## Mismatch combinations:
    ## AABB, AABC, ABCD, 2*AABR_AB, 2*AARB_BA, 2*ARRB_BA, 2*ABCR_ABC, 2*ABRC_CAB
    mismatch <- pp[c("AABB","AABC","ABCD","AABR_AB","AABR_AB","AARB_BA","AARB_BA",
                     "ARRB_BA","ARRB_BA","ABCR_ABC","ABCR_ABC","ABRC_CAB","ABRC_CAB")]
    ## Partial combinations:
    ## RARB, AAAB, ABAC, 2*AABR_BA, 2*AARB_AB, BAAR, ABRA, ABBR, BARB, 2*ARRB_AB
    ## ARBR_AB, ARBR_BA, AAAR, 2*ABCR_ACB, 2*ABCR_CAB, 2*ABRC_ABC, 2*ABRC_ACB
    partial <- pp[c("RARB","AAAB","ABAC","AABR_BA","AABR_BA","AARB_AB","AARB_AB","BAAR",
                    "ABRA","ABBR","BARB","ARRB_AB","ARRB_AB","ARBR_AB","ARBR_BA","AAAR",
                    "ABCR_ACB","ABCR_ACB","ABCR_CAB","ABCR_CAB","ABRC_ABC","ABRC_ABC",
                    "ABRC_ACB","ABRC_ACB")]
    ## Match combinations:
    ## AARR, ABRR, ARAR, ARRA, ARRR, RARA, RARR, RRRR, AAAA, ABAB, ABAR, BARA, BABR, ABRB
    match <- pp[c("AARR","ABRR","ARAR","ARRA","ARRR","RARA","RARR",
                  "RRRR","AAAA","ABAB","ABAR","BARA","BABR","ABRB")]
    c(sum(mismatch),sum(partial),sum(match))
}


## Makes all permuations of a vector and returns a matrix // Based on the 'multicool' package of Prof. James M. Curran
permAll <- function(x){
  if(length(x)==1) return(x)
  xx = initMC(x)
  allPerm(xx)
}


## Computes and returns the expected value of the cell counts
dbExpect <- function(probs,theta=0,k=c(0,0,1),n=1,r=0,R=0,round=FALSE,na=TRUE,vector=FALSE,collapse=FALSE,wildcard=FALSE,no.wildcard=NULL,rare.allele=FALSE,no.rare.allele=NULL){
    if(length(theta)>1)
        return(lapply(theta,function(t)
                      dbExpect(probs=probs,theta=t,k=k,n=n,r=r,R=R,round=round,na=na,vector=vector,collapse=collapse,
                               wildcard=wildcard,no.wildcard=no.wildcard,rare.allele=rare.allele,no.rare.allele=no.rare.allele)))
    if(rare.allele & wildcard) stop("Only one of 'wildcard' and 'rare.allele' can be TRUE.\nOtherwise set 'no.wildcard' and 'no.rare.allele'.")
    if(wildcard) f <- "F"
    else if(rare.allele){
        if(any(c(unlist(r),unlist(R))>0)) f <- "R" ## If one or more thresholds or tail probabilities are set
        else f <- "" ## If no probability is assigned the rare alleles
    }
    else f <- ""
    if(!is.list(probs) && is.vector(probs)) return(get(paste(f,"Ps",sep=""))(probs,t=theta,k=k,r=r,R=R)) ## Previous line handles if more thetas are provided 
    ## probs is a list of vectors with each vector being
    ## the allele probabilities for a given locus
    S <- length(probs)
    if(any(unlist(r)>0)){
        if(!(length(r)==1 | length(r)==S)) stop("Length of 'r' must be equal to the number of loci.")
        if(length(r)==1) r <- rep(r,S)
    }
    if(any(unlist(R)>0)){
        if(!(length(unlist(R)) %in% c(1,2,S,2*S))) stop("Length of 'R' not correct.")
        if(length(R)==1 && is.list(R)) R <- replicate(S,R,simplify=FALSE)
        else if(!is.list(R)) R <- replicate(S,R,simplify=FALSE)
    }
    if(is.character(k)){
        if(!is.na(match(toupper(k),c("UN","AV","FS","FC","PC"))))
            k <- list("UN"=c(0,0,1),"FC"=c(0,1,3)/4,"AV"=c(0,1,1)/2,"PC"=c(0,1,0),"FS"=c(1,2,1)/4)[[toupper(k)]]
        else stop(paste(paste("The value of 'k' (",k,") is not defined - has to be vector c(k2,k1,k0) or string.",sep=""),
                        "Options are: 'UN' (unrelated), 'FS' (full-siblings), 'AV' (avuncular), 'FC' (first cousins) or 'PC' (parent-child)",sep="\n"))
    }
    if(length(r)>1){
        if(length(R)!=length(r)) R <- replicate(length(r),R,simplify=FALSE)
        p <- lapply(as.list(1:S),function(z) RPs(probs[[z]],r=r[z],t=theta,k=k,R=R[[z]]))
    }
    else p <- lapply(probs,get(paste(f,"Ps",sep="")),t=theta,k=k,r=r,R=R)
    q <- do.call("cbind",p)
    res <- rare(q)
    ## In case of mixed matching
    if(!is.null(no.wildcard) | !is.null(no.rare.allele)){
        if(is.null(no.wildcard)) no.wildcard <- 0
        if(is.null(no.rare.allele)) no.rare.allele <- 0
        if((no.wildcard+no.rare.allele)>S) stop("Too many wildcard and rare allele loci compared to total number of loci")
        if(no.wildcard==0 & no.rare.allele==0) return(dbExpect(probs=probs,theta=theta,k=k,n=n,r=r,R=R,round=round,na=na,vector=vector,collapse=collapse,wildcard=FALSE,rare.allele=FALSE))
        else if(no.wildcard==S || no.rare.allele==S){ ## If w+r equals L // No exact matching
            if(no.wildcard==S) ## If all wildcard -> no rare
                return(dbExpect(probs=probs,theta=theta,k=k,n=n,r=r,R=R,round=round,na=na,vector=vector,collapse=collapse,wildcard=TRUE,rare.allele=FALSE))
            else ## 
                return(dbExpect(probs=probs,theta=theta,k=k,n=n,r=r,R=R,round=round,na=na,vector=vector,collapse=collapse,wildcard=FALSE,rare.allele=TRUE))
        }
        else{ ## If mixture 
            M <- permAll(rep(0:2,c(no.wildcard,no.rare.allele,S-(no.wildcard+no.rare.allele))))
            nM <- nrow(M)
            PPs <- do.call("cbind",lapply(1:length(probs),function(x) Ps(p=probs[[x]],t=theta,k=k,r=r[x],R=R[[x]])))
            FFs <- do.call("cbind",lapply(1:length(probs),function(x) FPs(p=probs[[x]],t=theta,k=k,r=r[x],R=R[[x]])))
            RRs <- do.call("cbind",lapply(1:length(probs),function(x) RPs(p=probs[[x]],t=theta,k=k,r=r[x],R=R[[x]])))
            res <- rare(cbind(FFs[,M[1,]==0],RRs[,M[1,]==1],PPs[,M[1,]==2]))/nM
            for(i in 2:nM) res <- res + rare(cbind(FFs[,M[i,]==0],RRs[,M[i,]==1],PPs[,M[i,]==2]))/nM
        }
    }
    if(n>1) N <- choose(n,2) else N <- 1
    if(round) res <- round(res*N)
    else res <- res*N
    if(na) res[!up.tri(res)] <- NA
    if(collapse) return(dbCollapse(res))
    if(vector){
        res <- t(res)[up.tri(res)]
        names(res) <- dbCats(S,vector=TRUE)
    }
    else dimnames(res) <- list(match=0:S,partial=0:S)
    res
}

