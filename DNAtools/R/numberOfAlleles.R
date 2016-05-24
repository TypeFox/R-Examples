## counts how many copies there is of each alpha-vector and returns the
## unique vectors together with their counts
permSummary.noa <- function(x){
  y <- table(sapply(x,num2pow))
  z <- lapply(strsplit(names(y),""),pow2num)
  list(type=z,terms=as.numeric(y))
}

## computes sum[i,...,j]^{different} P(A[i]...A[j])
Sa.noa <- function(a,p){
  if(length(a)==1) return(sum((p^a)))
  else{
    perm <- permSummary.noa(updateAlpha(a))
    return(Sa.noa(a[length(a)],p=p)*Sa.noa(a[-length(a)],p=p) - sum(perm$terms*unlist(lapply(perm$type,Sa.noa,p=p))))
  }
}

## computes the sum over sum[] P(AAAAA) including theta correction
Sab.noa <- function(a=rep(0,length(b)),b,t,p){
  if(all(t==0)) return(Sa.noa(a+b,p=p))
  if(length(b)==0){
    ## cat(paste("Sa(a=c(",paste(a,collapse=","),"),p=p)*(1/(1-t))\n\n",sep=""))
    return(Sa.noa(a,p)*(1/(1-t)))
  }
  else if(any(b==0)) return(Sab.noa(a,b[b!=0],t,p))
  else{
    if(b[length(b)]==1){
      return((1-t)*Sab.noa(a+ek(a,length(b)),b[-length(b)],t,p))
    }
    else{
      return((b[length(b)]-1)*t*Sab.noa(a,b-ek(b,length(b)),t,p) + (1-t)*Sab.noa(a+ek(a,length(b)),b-ek(b,length(b)),t,p))
    }
  }
}

## Converts char power representation to numeric vector
pow2num <- function(x){
  x <- sort(unlist(strsplit(paste(x),"")))
  x[grepl("[a-z]",x)] <- which(is.element(letters,x[grepl("[a-z]",x)]))+9 ## string to numeric; a=10, b=11, c=12, ...
  as.numeric(x)
}

## Converts numeric to char power
num2pow <- function(x){
  x <- x[x!=0]
  x[x>9] <- letters[x[x>9]-9]
  paste(sort(x),collapse="")
}

## If some of the derived functions for probability computations yield negative probabilities,
## then this is most likely due to numerical issues (arithmetics with very small numbers)
## This functions collapse the 'min' rarest alleles into one 'hyper allele':
collapseSmallest <- function(x,min=2){  
  x <- sort(x)
  smx <- sum(x[1:min])
  c(small=smx,x[-(1:min)])
}

## P[N(m)=n]
PNm.locus <- function(m=2,freq,theta=0,noa.tabs=NULL,maxA=NULL){
  if(2*m>length(freq)) maxA <- length(freq) ## Limited number of different alleles
  ## Denominator i P(A|AAAAAA) which is 1+(n-1)*theta. For 2m alleles the maximum P(A|..AA..) is 2m-1 alleles,
  ## hence 1+[1...2m-2]*theta which is 1+[1...2*(m-1)]*theta 
  denominator <- ifelse(m==1,1,prod(1+(1:(2*(m-1)))*theta))
  ## Create data frame with relevant information
  if(is.null(noa.tabs)) noa.tabs <- noaTabs(m=m)
  noatabs <- cbind(data.frame(count=noa.tabs),
                   structure(as.data.frame(do.call("rbind",strsplit(names(noa.tabs),"\\."))),.Names=c("noa","power")))
  if(!is.null(maxA)) noatabs <- noatabs[as.numeric(paste(noatabs$noa))<=maxA,] ## R CMD check failed this line: subset(noatabs,as.numeric(paste(noa))<=maxA)
  noatabs$p <- apply(noatabs,1,function(y) as.numeric(y[1])*Sab.noa(b=pow2num(y[3]),p=freq,t=theta))/denominator
  ## Output the result
  res <- aggregate(p~noa,data=noatabs,sum)
  ## If any probabilities are negative. This is most likely caused by very small allele frequencies for a number of alleles.
  ## In order to avoid this, the smallest alleles are collapsed into a 'hyper allele' covering the 'rare alleles':
  if(any(res$p<0)){
    warning("Negative probabilities: Collapsing the two smallest probabilities into a 'hyper allele'")
    return(PNm.locus(m=m,freq=collapseSmallest(freq),theta=theta,noa.tabs=noa.tabs,maxA=maxA))
  }
  res$noa <- as.numeric(paste(res$noa))
  if(!is.null(maxA)) res <- rbind(res,data.frame("noa"=(maxA+1):(2*m),p=0))
  res <- res[order(res$noa),]
  rownames(res) <- 1:nrow(res)
  res
}

## Recursion over loci in order to obtain result for P(N(m)=n) over all loci
PNL <- function(x){
  if(is.data.frame(x)) return(structure(x$p,.Names=x$noa))
  if(!is.list(x)) stop("Input must either be a data frame or a list of data frames")
  L <- length(x)
  noa <- nrow(x[[1]])
  m <- noa/2
  P <- replicate(L,rep(0,L*noa),simplify=FALSE)
  P[[1]][1:noa] <- x[[1]]$p
  for(l in 2:L){
    for(nn in l:(noa*l)){
      I <- min(noa,nn-1)
      P[[l]][nn] <- sum(P[[l-1]][nn-(1:I)]*x[[l]]$p[1:I])
    }
  }
  names(P[[L]]) <- 1:length(P[[L]])
  P[[L]]
}

## Wrapper function

pNoA <- function(probs, m=2, theta = 0, noa.tabs=NULL,locuswise=FALSE){
  ### probs: List of vectors with allele frequencies
  ### m: The number of contributors
  ### noa.tabs: alpha_m-vectors and their weights, c(alpha_m). If computed previously this may speed up computation
  ## If no noa.tabs are provided, compute these using noaTabs-function
  if(is.null(noa.tabs)) noa.tabs <- noaTabs(m=m)
  PNloci <- lapply(probs,PNm.locus,m=m,theta=theta,noa.tabs=noa.tabs)
  if(locuswise) return(lapply(PNloci,function(x) structure(x$p,.Names=1:nrow(x))))
  PNL(PNloci)
}

## Another name for same function
p.numberofalleles <- pNoA

### Functions for assessing the drop-out probability
PNdrop <- function(pD,probs,n0,m,loci){
  i <- 0:((2*m*loci)-n0)
  ((1-pD)^n0)*sum(choose(n0+i,i)*probs[n0+i]*(pD^i))
}


### Input the number of observed alleles (n0), number of contributors (m),
### and a vector of probabilities OR a list of probability vectors
estimatePD <- function(n0,m,pnoa=NULL,probs=NULL,theta=0,noa.tabs=NULL,locuswise=FALSE){
  if(locuswise){
    pDs <- structure(rep(NA,length(probs)),.Names=names(probs))
    pDs[n0==0] <- 1
  }
  else if(n0==0) return(c("P(D)"=1))
  if(is.null(pnoa)){
    if(is.null(probs)) stop("Either 'pNoA' or 'probs' must be provided")
    else{
      if(is.null(noa.tabs)) noa.tabs <- noaTabs(m=m)
      pnoa <- pNoA(probs,m=m,theta=theta, noa.tabs=noa.tabs, locuswise=locuswise)
    }
  }
  if(locuswise & any(n0>2*m))
    stop(paste("One of the number of observed alleles (",(n0[n0>2*m])[1],") is bigger than the maximum for a ",m,"-person mixture (",2*m,").",sep=""))
  else if(!locuswise & any(n0>length(pnoa)))
    stop(paste("The number of observed alleles (",n0,") is bigger than the maximum for a ",m,"-person mixture (",length(pnoa),").",sep=""))
  if(locuswise){
    if(any(n0!=0)) pDs[n0!=0] <- mapply(function(x,y)optimise(PNdrop,m=m,probs=x,n0=y,loci=1,lower=0,upper=1,maximum=TRUE)$maximum,x=pnoa[n0!=0],y=n0[n0!=0])
    return(pDs)
  }
  else{
    loci <- length(pnoa)/(2*m)
    return(structure(optimise(PNdrop,m=m,probs=pnoa,n0=n0,loci=loci,lower=0,upper=1,maximum=TRUE)$maximum,.Names=c("p(D)")))
  }
}

### Workhorse of the package - compute the alpha-vectors by recursion

locusrecur <- function(alpha,weight=1,m){
  if(sum(alpha)>2*m) stop("There is more terms and 2*m with m being the number of profiles")
  if(sum(alpha)==2*m) return(c(num2pow(alpha),weight)) ## Return alpha and weight
  if(is.null(alpha)) return(mapply(locusrecur,alpha=list(c(1,1),c(2)),weight=rep(1,2),m=m,SIMPLIFY=FALSE))
  la <- length(alpha)
  ## Addition of het genotype
  ## First profile is hom - hence can't share two alleles with het
  if(la>=2){ hets2 <- cbind(permAll(c(rep(0,la-2),1,1)),0,0); nh2 <- nrow(hets2) }  ## weight: 2
  else{ hets2 <- NULL; nh2 <- 0 }
  ## Addition of one allele ok
  hets1 <- cbind(permAll(c(rep(0,la-1),1)),1,0) ## weight: 2
  hets0 <- matrix(c(rep(0,la),1,1),nrow=1) ## weight: 1
  hets <- rbind(hets2,hets1,hets0)
  weights <- rep(c(2,2,1),c(nh2,nrow(hets1),nrow(hets0)))
  ## Addition of hom genotype
  homs <- permAll(c(rep(0,la),2))
  weights <- weight*c(weights,rep(1,nrow(homs)))
  ## hets <- permAll(c(rep(0,length(alpha)),1,1)) ## Addition of het genotype
  ## weights <- 1+apply(hets,1,function(z) all(z[length(z)-(1:0)]==0)) ## fully overlapping hets are counted twice
  ## homs <- permAll(c(rep(0,length(alpha)),2)) ## Addition of hom genotype; weight: 1
  ## weights <- weight*c(weights,rep(1,nrow(homs)))
  alphahets <- t(t(hets)+c(alpha,0,0))
  alphahoms <- t(t(homs)+c(alpha,0))
  alphas <- c(split(alphahets,1:nrow(alphahets)),split(alphahoms,1:nrow(alphahoms)))
  alphas <- unlist(lapply(alphas,num2pow))
  aws <- aggregate(w~a,data=data.frame(w=weights,a=alphas),sum)
  alphas <- sapply(aws$a,pow2num,simplify=FALSE)
  if(any(unlist(lapply(alphas,sum))>2*(m))) browser()
  mapply(locusrecur,alpha=alphas,weight=aws$w,m=m,SIMPLIFY=FALSE)
}

noaTabs <- function(alpha=NULL,m=2,weight=1){
  if(m==1) return(as.array(c("1.2"=1,"2.11"=1),dim=2))
  if(is.array(alpha)){
    alphas <- lapply(strsplit(names(alpha),"\\."),function(x) pow2num(x[2]))
    weights <- as.list(as.numeric(alpha))
    rL <- unlist(mapply(noaTabs,alpha=alphas,weight=weights,m=m,SIMPLIFY=FALSE))
    rL <- aggregate(count~type,data=data.frame(type=names(rL),count=as.numeric(rL)),sum)
    rL <- structure(rL$count,.Dim=nrow(rL),.Names=paste(rL$type))
    return(rL[sort.list(nchar(names(rL)))])
  }
  rL <- locusrecur(alpha=alpha,m=m,weight=weight)
  rL <- structure(as.data.frame(matrix(unlist(rL),ncol=2,byrow=TRUE)),.Names=c("type","count"))
  rL <- aggregate(count~type,data=within(rL,{count <- as.numeric(paste(count)); type <- paste(type)}),sum)
  rL <- structure(rL$count,.Dim=nrow(rL),.Names=rL$type)
  names(rL) <- paste(nchar(names(rL)),names(rL),sep=".")
  rL <- rL[sort.list(nchar(names(rL)))]
  rL
}

### Code for evaluating P(m|n0)

pContrib.locus <- function(prob=NULL, m.prior=NULL, m.max=8, pnoa.locus=NULL, theta=0){
  if(is.null(m.prior)) m.prior <- rep(1/m.max,m.max)
  if(length(m.prior)!=m.max) m.max <- length(m.prior)
  if(!is.null(pnoa.locus)){
    if(!is.list(pnoa.locus)) stop("'pnoa.locus' must be a named list of locus probabilities, e.g. list('1'=locus.m1,'2'=locus.m2)")
    if(!all(1:m.max %in% names(pnoa.locus)) & is.null(prob)) stop("If not all locuswise probabilities are supplied for 1,...,m.max, then a vector of allele frequencies must be supplied.") 
    PNOA.locus <- structure(replicate(m.max,NULL,simplify=FALSE),.Names=1:m.max)
    for(i in 1:m.max){
      if(i%in%names(pnoa.locus)) PNOA.locus[[i]] <- c(pnoa.locus[[paste(i)]],rep(0,2*(m.max-i)))
      else PNOA.locus[[i]] <- c(PNm.locus(m=i,freq=prob,theta=theta)$p,rep(0,2*(m.max-i)))
    }
  }
  else PNOA.locus <- lapply(1:m.max,function(i) c(PNm.locus(m=i,prob,theta=theta)$p,rep(0,2*(m.max-i))))
  PNOA.locus <- do.call("rbind",PNOA.locus)
  ## Compute the posterior: P(m|n) = P(n|m)*P(m)/P(n); where the computations are done for n=1,...,2*m.max
  ## The denominator in vector format: P(n) = (P(n=1),..,P(n=2*m.max)) = (sum{j=1}^m.max P(n=1|m=j)*P(m=j), ..., sum{j=1}^m.max P(n=2*m.max|m=j)*P(m=j))
  denominator <- as.numeric(m.prior%*%PNOA.locus)
  ## Compute numerator: a vector of vectors (matrix): [P(n=1|m=1)*p(m=1) ... P(n=1|m=M)*P(m=M); .... ; P(n=2M|m=1)*p(m=1) ... P(n=2M|m=M)*P(m=M)]
  numerator <- t(PNOA.locus)*m.prior
  ## Return a matrix of posterior probabilities
  structure(numerator/denominator,.Dimnames=list("n0"=1:(2*m.max),"P(m|n0)\n   m"=1:m.max))
}

### Computes the posterior probabilities of m given a specific vector of observed locus counts
pContrib <- function(n0, probs=NULL, m.prior=rep(1/m.max,m.max), m.max=8, pnoa=NULL, theta=0){
  if(length(m.prior)!=m.max) m.max <- length(m.prior)
  if(!is.null(pnoa)){
    if(!is.list(pnoa)) stop("'pNoA' must be a named list of locus probabilities, e.g. list('1'=list(locus1.m1, locus2.m1, ...),'2'=list(locus1.m2, ...)")
    if(!all(1:m.max %in% names(pnoa)) & is.null(probs)) stop("If not all locuswise probabilities are supplied for 1,...,m.max, then a list of vectors with allele frequencies must be supplied.")    
    PNOA <- structure(replicate(m.max,NULL,simplify=FALSE),.Names=1:m.max)
    for(i in 1:m.max){
      if(i%in%names(pnoa)) PNOA[[i]] <- lapply(pnoa[[paste(i)]],function(x,ii) c(x,rep(0,2*(m.max-ii))),ii=i)
      else PNOA[[i]] <- lapply(suppressWarnings(pNoA(probs=probs, m=i, theta=theta, locuswise=TRUE)),function(x,ii) c(x,rep(0,2*(m.max-ii))),ii=i)
    }
  }
  else PNOA <- lapply(1:m.max,function(i) lapply(suppressWarnings(pNoA(probs=probs,m=i,theta=theta,locuswise=TRUE)),function(x,ii) c(x,rep(0,2*(m.max-ii))),ii=i))
  ## Extract the relevant entries from each locus: P(n_l|m=1),...,P(n_l|m=M).
  ## This is done by forming [P(n_1|m=1) P(n_2|m=1) ... P(n_L|m=1); ... ; P(n_1|m=M) ... P(n_L|m=M)]
  ## Furthermore, the numerator is evaluated by P(m)*prod_{l=1}^L P(n_0,l|m), where the product is evaluated from extracted entries, cf. above.
  pn0.m.prod <- unlist(lapply(PNOA,function(pnm) prod(mapply(function(pn,n)pn[n],pn=pnm,n=n0))))
  ## The numerator is then [P(m=1)*prod_{l=1}^L P(n_0,l|m=1) ... P(m=M)*prod_{l=1}^L P(n_0,l|m=M)]:
  numerator <- pn0.m.prod*m.prior ## holds P(m|n0)*P(m) for m=1,...,length(m.prior)
  ## The denominator is simply the sum  over m in the numerator:
  structure(numerator/sum(numerator),.Names=paste(paste("P(m=",1:m.max,sep=""),"|n0)",sep=""))
}
